import argparse
import uproot
import subprocess
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nHessians", type=int, required=True)
parser.add_argument("-p", "--pdf", type=str, required=False)
parser.add_argument("-v", "--variable", type=str, required=False)
args = parser.parse_args()

def runVar(idx, var, shift):
    new_text= []
    label = f"pdf{idx}{shift}"
    print("label is", label)
    with open("WGen.txt") as card:
        text = card.readlines()
        for line in text:
            if "kmax" in line:
                new_text.append("kmax 2\n")
            elif "shapes" in line:
                if "data_obs" in line:
                    new_text += line
                else:
                    temp =  line.replace(f"{var}_mp", f"{var}_{label}_mp")
                    new_text += temp.replace(f"{var}_mn", f"{var}_{label}_mn")
            elif "group" in line or ("shape" in line and not any([x in line for x in ("massShift100MeV", "lumi2016_13TeV")])):
                continue
            else:
                new_text += line
    cardname = f"WGen_{label}.txt"
    outfile = f"fitresults_{label}.root"
    with open(cardname, "w") as newcard:
        newcard.write("".join(new_text))
    subprocess.call(["text2hdf5.py", "--X-allow-no-signal", cardname])
    subprocess.call(["combinetf.py", cardname.replace("txt", "hdf5"), f"-o {outfile}"])
    return outfile

def getMassShift(rf):
    return rf["fitresults"]["massShift100MeV"].array(library="np")[0]

upvals = []
downvals = []
for i in range(1, args.nHessians+1):
    outu = runVar(i, args.variable, "Up")
    rf = uproot.open(outu)
    up = getMassShift(rf)

    outd = runVar(i, args.variable, "Down")
    rf = uproot.open(outd)
    down = getMassShift(rf)

    upvals.append(max(up, down))
    downvals.append(min(up, down))

print(upvals)
print(downvals)
scale = 1 if args.pdf != "ct18" else 0.6079
print("Scaling uncertainty by", scale)
upunc = np.sqrt(np.sum([x**2 for x in filter(lambda x: x > 0, upvals)]))*scale
downunc = np.sqrt(np.sum([x**2 for x in filter(lambda x: x < 0, downvals)]))*scale
print("Uncertinaty: +%0.2f/-%02f" % (upunc, downunc))
