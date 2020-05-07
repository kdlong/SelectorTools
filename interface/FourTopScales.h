#ifndef _FOURTOPSCALES_H_
#define _FOURTOPSCALES_H_

float getTruePUw(int nvtx) {
    if (nvtx >= 0.000000 && nvtx < 1.000000) return 0.366077;
    if (nvtx >= 1.000000 && nvtx < 2.000000) return 0.893925;
    if (nvtx >= 2.000000 && nvtx < 3.000000) return 1.197716;
    if (nvtx >= 3.000000 && nvtx < 4.000000) return 0.962699;
    if (nvtx >= 4.000000 && nvtx < 5.000000) return 1.120976;
    if (nvtx >= 5.000000 && nvtx < 6.000000) return 1.164859;
    if (nvtx >= 6.000000 && nvtx < 7.000000) return 0.795599;
    if (nvtx >= 7.000000 && nvtx < 8.000000) return 0.495824;
    if (nvtx >= 8.000000 && nvtx < 9.000000) return 0.742182;
    if (nvtx >= 9.000000 && nvtx < 10.000000) return 0.878856;
    if (nvtx >= 10.000000 && nvtx < 11.000000) return 0.964232;
    if (nvtx >= 11.000000 && nvtx < 12.000000) return 1.072499;
    if (nvtx >= 12.000000 && nvtx < 13.000000) return 1.125335;
    if (nvtx >= 13.000000 && nvtx < 14.000000) return 1.176027;
    if (nvtx >= 14.000000 && nvtx < 15.000000) return 1.202083;
    if (nvtx >= 15.000000 && nvtx < 16.000000) return 1.207643;
    if (nvtx >= 16.000000 && nvtx < 17.000000) return 1.200176;
    if (nvtx >= 17.000000 && nvtx < 18.000000) return 1.182682;
    if (nvtx >= 18.000000 && nvtx < 19.000000) return 1.143998;
    if (nvtx >= 19.000000 && nvtx < 20.000000) return 1.096632;
    if (nvtx >= 20.000000 && nvtx < 21.000000) return 1.065602;
    if (nvtx >= 21.000000 && nvtx < 22.000000) return 1.051166;
    if (nvtx >= 22.000000 && nvtx < 23.000000) return 1.051600;
    if (nvtx >= 23.000000 && nvtx < 24.000000) return 1.050630;
    if (nvtx >= 24.000000 && nvtx < 25.000000) return 1.049862;
    if (nvtx >= 25.000000 && nvtx < 26.000000) return 1.058173;
    if (nvtx >= 26.000000 && nvtx < 27.000000) return 1.072155;
    if (nvtx >= 27.000000 && nvtx < 28.000000) return 1.083030;
    if (nvtx >= 28.000000 && nvtx < 29.000000) return 1.095693;
    if (nvtx >= 29.000000 && nvtx < 30.000000) return 1.107871;
    if (nvtx >= 30.000000 && nvtx < 31.000000) return 1.094621;
    if (nvtx >= 31.000000 && nvtx < 32.000000) return 1.082620;
    if (nvtx >= 32.000000 && nvtx < 33.000000) return 1.041247;
    if (nvtx >= 33.000000 && nvtx < 34.000000) return 0.985752;
    if (nvtx >= 34.000000 && nvtx < 35.000000) return 0.910807;
    if (nvtx >= 35.000000 && nvtx < 36.000000) return 0.820923;
    if (nvtx >= 36.000000 && nvtx < 37.000000) return 0.716787;
    if (nvtx >= 37.000000 && nvtx < 38.000000) return 0.610013;
    if (nvtx >= 38.000000 && nvtx < 39.000000) return 0.503118;
    if (nvtx >= 39.000000 && nvtx < 40.000000) return 0.404841;
    if (nvtx >= 40.000000 && nvtx < 41.000000) return 0.309195;
    if (nvtx >= 41.000000 && nvtx < 42.000000) return 0.227920;
    if (nvtx >= 42.000000 && nvtx < 43.000000) return 0.163690;
    if (nvtx >= 43.000000 && nvtx < 44.000000) return 0.113180;
    if (nvtx >= 44.000000 && nvtx < 45.000000) return 0.077300;
    if (nvtx >= 45.000000 && nvtx < 46.000000) return 0.050922;
    if (nvtx >= 46.000000 && nvtx < 47.000000) return 0.031894;
    if (nvtx >= 47.000000 && nvtx < 48.000000) return 0.020094;
    if (nvtx >= 48.000000 && nvtx < 49.000000) return 0.012263;
    if (nvtx >= 49.000000 && nvtx < 50.000000) return 0.007426;
    if (nvtx >= 50.000000 && nvtx < 51.000000) return 0.004380;
    if (nvtx >= 51.000000 && nvtx < 52.000000) return 0.002608;
    if (nvtx >= 52.000000 && nvtx < 53.000000) return 0.001566;
    if (nvtx >= 53.000000 && nvtx < 54.000000) return 0.000971;
    if (nvtx >= 54.000000 && nvtx < 55.000000) return 0.000729;
    if (nvtx >= 55.000000 && nvtx < 56.000000) return 0.000673;
    if (nvtx >= 56.000000 && nvtx < 57.000000) return 0.000730;
    if (nvtx >= 57.000000 && nvtx < 58.000000) return 0.000949;
    if (nvtx >= 58.000000 && nvtx < 59.000000) return 0.001355;
    if (nvtx >= 59.000000 && nvtx < 60.000000) return 0.001894;
    if (nvtx >= 60.000000 && nvtx < 61.000000) return 0.003082;
    if (nvtx >= 61.000000 && nvtx < 62.000000) return 0.004097;
    if (nvtx >= 62.000000 && nvtx < 63.000000) return 0.004874;
    if (nvtx >= 63.000000 && nvtx < 64.000000) return 0.005256;
    if (nvtx >= 64.000000 && nvtx < 65.000000) return 0.005785;
    if (nvtx >= 65.000000 && nvtx < 66.000000) return 0.005515;
    if (nvtx >= 66.000000 && nvtx < 67.000000) return 0.005000;
    if (nvtx >= 67.000000 && nvtx < 68.000000) return 0.004410;
    if (nvtx >= 68.000000 && nvtx < 69.000000) return 0.004012;
    if (nvtx >= 69.000000 && nvtx < 70.000000) return 0.003548;
    if (nvtx >= 70.000000 && nvtx < 71.000000) return 0.003108;
    if (nvtx >= 71.000000 && nvtx < 72.000000) return 0.002702;
    if (nvtx >= 72.000000 && nvtx < 73.000000) return 0.002337;
    if (nvtx >= 73.000000 && nvtx < 74.000000) return 0.002025;
    if (nvtx >= 74.000000 && nvtx < 75.000000) return 0.001723;
    if (nvtx >= 75.000000 && nvtx < 76.000000) return 0.000000;
    if (nvtx >= 76.000000 && nvtx < 77.000000) return 0.000000;
    if (nvtx >= 77.000000 && nvtx < 78.000000) return 0.000000;
    if (nvtx >= 78.000000 && nvtx < 79.000000) return 0.000000;
    if (nvtx >= 79.000000 && nvtx < 80.000000) return 0.000000;
    return 0.;
}

float getTruePUw_Moriond(int nvtx) {
    // below derived using SS baby
    // if (nvtx>=0.000000 && nvtx<1.000000) return 0.339493;
    // if (nvtx>=1.000000 && nvtx<2.000000) return 1.587625;
    // if (nvtx>=2.000000 && nvtx<3.000000) return 1.193400;
    // if (nvtx>=3.000000 && nvtx<4.000000) return 0.888480;
    // if (nvtx>=4.000000 && nvtx<5.000000) return 1.106730;
    // if (nvtx>=5.000000 && nvtx<6.000000) return 0.947135;
    // if (nvtx>=6.000000 && nvtx<7.000000) return 0.674775;
    // if (nvtx>=7.000000 && nvtx<8.000000) return 0.433799;
    // if (nvtx>=8.000000 && nvtx<9.000000) return 0.609507;
    // if (nvtx>=9.000000 && nvtx<10.000000) return 0.765006;
    // if (nvtx>=10.000000 && nvtx<11.000000) return 0.851999;
    // if (nvtx>=11.000000 && nvtx<12.000000) return 0.943593;
    // if (nvtx>=12.000000 && nvtx<13.000000) return 1.022762;
    // if (nvtx>=13.000000 && nvtx<14.000000) return 1.050225;
    // if (nvtx>=14.000000 && nvtx<15.000000) return 1.107856;
    // if (nvtx>=15.000000 && nvtx<16.000000) return 1.104652;
    // if (nvtx>=16.000000 && nvtx<17.000000) return 1.090875;
    // if (nvtx>=17.000000 && nvtx<18.000000) return 1.116378;
    // if (nvtx>=18.000000 && nvtx<19.000000) return 1.056206;
    // if (nvtx>=19.000000 && nvtx<20.000000) return 1.041334;
    // if (nvtx>=20.000000 && nvtx<21.000000) return 1.025138;
    // if (nvtx>=21.000000 && nvtx<22.000000) return 1.016013;
    // if (nvtx>=22.000000 && nvtx<23.000000) return 1.019632;
    // if (nvtx>=23.000000 && nvtx<24.000000) return 1.043251;
    // if (nvtx>=24.000000 && nvtx<25.000000) return 1.015568;
    // if (nvtx>=25.000000 && nvtx<26.000000) return 1.070650;
    // if (nvtx>=26.000000 && nvtx<27.000000) return 1.093817;
    // if (nvtx>=27.000000 && nvtx<28.000000) return 1.084384;
    // if (nvtx>=28.000000 && nvtx<29.000000) return 1.126802;
    // if (nvtx>=29.000000 && nvtx<30.000000) return 1.154778;
    // if (nvtx>=30.000000 && nvtx<31.000000) return 1.143681;
    // if (nvtx>=31.000000 && nvtx<32.000000) return 1.125193;
    // if (nvtx>=32.000000 && nvtx<33.000000) return 1.095793;
    // if (nvtx>=33.000000 && nvtx<34.000000) return 1.047938;
    // if (nvtx>=34.000000 && nvtx<35.000000) return 1.019224;
    // if (nvtx>=35.000000 && nvtx<36.000000) return 0.900380;
    // if (nvtx>=36.000000 && nvtx<37.000000) return 0.774526;
    // if (nvtx>=37.000000 && nvtx<38.000000) return 0.688022;
    // if (nvtx>=38.000000 && nvtx<39.000000) return 0.550549;
    // if (nvtx>=39.000000 && nvtx<40.000000) return 0.472619;
    // if (nvtx>=40.000000 && nvtx<41.000000) return 0.360133;
    // if (nvtx>=41.000000 && nvtx<42.000000) return 0.271264;
    // if (nvtx>=42.000000 && nvtx<43.000000) return 0.192419;
    // if (nvtx>=43.000000 && nvtx<44.000000) return 0.136818;
    // if (nvtx>=44.000000 && nvtx<45.000000) return 0.093556;
    // if (nvtx>=45.000000 && nvtx<46.000000) return 0.065444;
    // if (nvtx>=46.000000 && nvtx<47.000000) return 0.041100;
    // if (nvtx>=47.000000 && nvtx<48.000000) return 0.027266;
    // if (nvtx>=48.000000 && nvtx<49.000000) return 0.015408;
    // if (nvtx>=49.000000 && nvtx<50.000000) return 0.009769;
    // below derived using CMS3
    if (nvtx >= 0.000000 && nvtx < 1.000000) return 0.456317;
    if (nvtx >= 1.000000 && nvtx < 2.000000) return 0.868822;
    if (nvtx >= 2.000000 && nvtx < 3.000000) return 1.060915;
    if (nvtx >= 3.000000 && nvtx < 4.000000) return 0.804616;
    if (nvtx >= 4.000000 && nvtx < 5.000000) return 1.279666;
    if (nvtx >= 5.000000 && nvtx < 6.000000) return 1.377812;
    if (nvtx >= 6.000000 && nvtx < 7.000000) return 0.704106;
    if (nvtx >= 7.000000 && nvtx < 8.000000) return 0.494254;
    if (nvtx >= 8.000000 && nvtx < 9.000000) return 0.718032;
    if (nvtx >= 9.000000 && nvtx < 10.000000) return 0.853469;
    if (nvtx >= 10.000000 && nvtx < 11.000000) return 0.938278;
    if (nvtx >= 11.000000 && nvtx < 12.000000) return 1.057068;
    if (nvtx >= 12.000000 && nvtx < 13.000000) return 1.122293;
    if (nvtx >= 13.000000 && nvtx < 14.000000) return 1.153207;
    if (nvtx >= 14.000000 && nvtx < 15.000000) return 1.205728;
    if (nvtx >= 15.000000 && nvtx < 16.000000) return 1.198719;
    if (nvtx >= 16.000000 && nvtx < 17.000000) return 1.184970;
    if (nvtx >= 17.000000 && nvtx < 18.000000) return 1.176247;
    if (nvtx >= 18.000000 && nvtx < 19.000000) return 1.127847;
    if (nvtx >= 19.000000 && nvtx < 20.000000) return 1.083409;
    if (nvtx >= 20.000000 && nvtx < 21.000000) return 1.063132;
    if (nvtx >= 21.000000 && nvtx < 22.000000) return 1.042620;
    if (nvtx >= 22.000000 && nvtx < 23.000000) return 1.045930;
    if (nvtx >= 23.000000 && nvtx < 24.000000) return 1.047821;
    if (nvtx >= 24.000000 && nvtx < 25.000000) return 1.042626;
    if (nvtx >= 25.000000 && nvtx < 26.000000) return 1.056079;
    if (nvtx >= 26.000000 && nvtx < 27.000000) return 1.068642;
    if (nvtx >= 27.000000 && nvtx < 28.000000) return 1.077145;
    if (nvtx >= 28.000000 && nvtx < 29.000000) return 1.090889;
    if (nvtx >= 29.000000 && nvtx < 30.000000) return 1.102057;
    if (nvtx >= 30.000000 && nvtx < 31.000000) return 1.089665;
    if (nvtx >= 31.000000 && nvtx < 32.000000) return 1.066530;
    if (nvtx >= 32.000000 && nvtx < 33.000000) return 1.037292;
    if (nvtx >= 33.000000 && nvtx < 34.000000) return 0.969761;
    if (nvtx >= 34.000000 && nvtx < 35.000000) return 0.900420;
    if (nvtx >= 35.000000 && nvtx < 36.000000) return 0.814996;
    if (nvtx >= 36.000000 && nvtx < 37.000000) return 0.701494;
    if (nvtx >= 37.000000 && nvtx < 38.000000) return 0.601403;
    if (nvtx >= 38.000000 && nvtx < 39.000000) return 0.492305;
    if (nvtx >= 39.000000 && nvtx < 40.000000) return 0.400363;
    if (nvtx >= 40.000000 && nvtx < 41.000000) return 0.308498;
    if (nvtx >= 41.000000 && nvtx < 42.000000) return 0.225170;
    if (nvtx >= 42.000000 && nvtx < 43.000000) return 0.162180;
    if (nvtx >= 43.000000 && nvtx < 44.000000) return 0.111863;
    if (nvtx >= 44.000000 && nvtx < 45.000000) return 0.077390;
    if (nvtx >= 45.000000 && nvtx < 46.000000) return 0.050697;
    if (nvtx >= 46.000000 && nvtx < 47.000000) return 0.031394;
    if (nvtx >= 47.000000 && nvtx < 48.000000) return 0.020052;
    if (nvtx >= 48.000000 && nvtx < 49.000000) return 0.012629;
    if (nvtx >= 49.000000 && nvtx < 50.000000) return 0.007430;
    return 0.;
}

float triggerScaleFactor(int pdgId1, int pdgId2, float pt1, float pt2,
                         float eta1, float eta2, float ht) {
    // return TotalTriggerSF(pdgId1, pt1, eta1, pdgId2, pt2, eta2, ht);
    // Using Matthieu's macro, so dummy 1 here
    return 1.0;  // FIXME

    if (ht > 300) {
        if ((abs(pdgId1) + abs(pdgId2)) == 22) return 1.;
        if ((abs(pdgId1) + abs(pdgId2)) == 26) return 0.985 * 0.985;
        if ((abs(pdgId1) + abs(pdgId2)) == 24) return 0.985;
    } else {
        if ((abs(pdgId1) + abs(pdgId2)) == 22) return 0.997 * 0.997 * 0.998;
        if ((abs(pdgId1) + abs(pdgId2)) == 26) return 0.982 * 0.985 * 0.973;
        if ((abs(pdgId1) + abs(pdgId2)) == 24) {
            if (abs(pdgId1) == 11) {
                if (pt1 > pt2)
                    return 0.997 * 0.985;
                else
                    return 0.997 * 0.982;
            } else {
                if (pt1 > pt2)
                    return 0.997 * 0.982;
                else
                    return 0.997 * 0.985;
            }
        }
    }
    return 0.;
}

float electronScaleFactorHighHT(float pt, float eta) {
    if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9158;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9820;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 1.0756;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9203;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 1.1124;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9177;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9499;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9710;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9370;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9904;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9210;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9472;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9927;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9443;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9785;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9213;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9515;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9830;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9480;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9628;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9212;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9483;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9845;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9480;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9483;
    if (pt >= 100 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9469;
    if (pt >= 100 && fabs(eta) >= 0.800 && fabs(eta) < 1.442) return 0.9429;
    if (pt >= 100 && fabs(eta) >= 1.442 && fabs(eta) < 1.566) return 0.9932;
    if (pt >= 100 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.9455;
    if (pt >= 100 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.9592;
    return 0.;
}
float electronScaleFactorLowHT(float pt, float eta) {
    if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9149;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9768;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 1.0781;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9169;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 1.1100;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9170;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9497;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9687;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9356;
    if (pt >= 20 && pt < 30 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9894;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9208;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9483;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9923;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9438;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9781;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9202;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9514;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9827;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9480;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9627;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 0.000 && fabs(eta) < 0.800)
        return 0.9207;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 0.800 && fabs(eta) < 1.442)
        return 0.9481;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 1.442 && fabs(eta) < 1.566)
        return 0.9848;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 1.566 && fabs(eta) < 2.000)
        return 0.9480;
    if (pt >= 50 && pt < 100 && fabs(eta) >= 2.000 && fabs(eta) < 2.500)
        return 0.9477;
    if (pt >= 100 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9472;
    if (pt >= 100 && fabs(eta) >= 0.800 && fabs(eta) < 1.442) return 0.9333;
    if (pt >= 100 && fabs(eta) >= 1.442 && fabs(eta) < 1.566) return 0.9934;
    if (pt >= 100 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.9383;
    if (pt >= 100 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.9597;
    return 0.;
}
float muonScaleFactor(float pt, float eta) {
    if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.900)
        return 0.9047;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 0.900 && fabs(eta) < 1.200)
        return 0.8860;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 1.200 && fabs(eta) < 2.100)
        return 0.8916;
    if (pt >= 10 && pt < 20 && fabs(eta) >= 2.100 && fabs(eta) < 2.400)
        return 0.8394;
    if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900)
        return 0.9430;
    if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200)
        return 0.9685;
    if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100)
        return 0.9741;
    if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400)
        return 0.8917;
    if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900)
        return 0.9707;
    if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200)
        return 0.9724;
    if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100)
        return 0.9777;
    if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400)
        return 0.9180;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900)
        return 0.9821;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200)
        return 0.9850;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100)
        return 0.9934;
    if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400)
        return 0.9389;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900)
        return 0.9854;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200)
        return 0.9861;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100)
        return 0.9968;
    if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400)
        return 0.9453;
    if (pt >= 50 && pt < 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900)
        return 0.9813;
    if (pt >= 50 && pt < 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200)
        return 0.9819;
    if (pt >= 50 && pt < 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100)
        return 0.9964;
    if (pt >= 50 && pt < 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400)
        return 0.9410;
    if (pt >= 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9830;
    if (pt >= 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9861;
    if (pt >= 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9994;
    if (pt >= 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9525;
    return 0.;
}

float electronGSF(float pt, float eta) {
    if (pt >= 20 && eta >= -2.500 && eta < -2.400) return 1.1703;
    if (pt >= 20 && eta >= -2.400 && eta < -2.300) return 1.0085;
    if (pt >= 20 && eta >= -2.300 && eta < -2.200) return 1.0105;
    if (pt >= 20 && eta >= -2.200 && eta < -2.000) return 1.0052;
    if (pt >= 20 && eta >= -2.000 && eta < -1.800) return 0.9979;
    if (pt >= 20 && eta >= -1.800 && eta < -1.630) return 0.9917;
    if (pt >= 20 && eta >= -1.630 && eta < -1.566) return 0.9865;
    if (pt >= 20 && eta >= -1.566 && eta < -1.444) return 0.9616;
    if (pt >= 20 && eta >= -1.444 && eta < -1.200) return 0.9867;
    if (pt >= 20 && eta >= -1.200 && eta < -1.000) return 0.9775;
    if (pt >= 20 && eta >= -1.000 && eta < -0.600) return 0.9694;
    if (pt >= 20 && eta >= -0.600 && eta < -0.400) return 0.9664;
    if (pt >= 20 && eta >= -0.400 && eta < -0.200) return 0.9633;
    if (pt >= 20 && eta >= -0.200 && eta < 0.000) return 0.9600;
    if (pt >= 20 && eta >= 0.000 && eta < 0.200) return 0.9662;
    if (pt >= 20 && eta >= 0.200 && eta < 0.400) return 0.9796;
    if (pt >= 20 && eta >= 0.400 && eta < 0.600) return 0.9766;
    if (pt >= 20 && eta >= 0.600 && eta < 1.000) return 0.9807;
    if (pt >= 20 && eta >= 1.000 && eta < 1.200) return 0.9867;
    if (pt >= 20 && eta >= 1.200 && eta < 1.444) return 0.9867;
    if (pt >= 20 && eta >= 1.444 && eta < 1.566) return 0.9707;
    if (pt >= 20 && eta >= 1.566 && eta < 1.630) return 0.9897;
    if (pt >= 20 && eta >= 1.630 && eta < 1.800) return 0.9959;
    if (pt >= 20 && eta >= 1.800 && eta < 2.000) return 0.9897;
    if (pt >= 20 && eta >= 2.000 && eta < 2.200) return 0.9949;
    if (pt >= 20 && eta >= 2.200 && eta < 2.300) return 0.9928;
    if (pt >= 20 && eta >= 2.300 && eta < 2.400) return 0.9666;
    if (pt >= 20 && eta >= 2.400 && eta < 2.500) return 0.8840;
    return 0.;
}

float trackingSF(float eta) {
    if (eta >= -2.400 && eta < -2.100) return 0.9879;
    if (eta >= -2.100 && eta < -1.600) return 0.9939;
    if (eta >= -1.600 && eta < -1.100) return 0.9970;
    if (eta >= -1.100 && eta < -0.600) return 0.9954;
    if (eta >= -0.600 && eta < 0.000) return 0.9937;
    if (eta >= 0.000 && eta < 0.600) return 0.9959;
    if (eta >= 0.600 && eta < 1.100) return 0.9976;
    if (eta >= 1.100 && eta < 1.600) return 0.9961;
    if (eta >= 1.600 && eta < 2.100) return 0.9930;
    if (eta >= 2.100 && eta < 2.400) return 0.9819;
    return 1.0;
}

float leptonScaleFactor(int pdgId, float pt, float eta, float ht) {
    if (abs(pdgId) == 13)
        return muonScaleFactor(pt, eta) * trackingSF(eta);
    else if (abs(pdgId) == 11) {
        if (ht > 300)
            return electronScaleFactorHighHT(pt, eta) * electronGSF(pt, eta);
        else
            return electronScaleFactorLowHT(pt, eta) * electronGSF(pt, eta);
    }
    return 0.;
}

float eventScaleFactor(int pdgId1, int pdgId2, float pt1, float pt2, float eta1,
                       float eta2, float ht) {
    return triggerScaleFactor(pdgId1, pdgId2, pt1, pt2, eta1, eta2, ht) *
           leptonScaleFactor(pdgId1, pt1, eta1, ht) *
           leptonScaleFactor(pdgId2, pt2, eta2, ht);
}

#endif
