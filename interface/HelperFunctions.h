#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <math.h>
#include "Analysis/VVAnalysis/interface/GoodParticle.h"
/* #include <pair> */

double JetSphericity(std::vector<GoodPart>& jets) {
    TMatrixDSym sphere(3);
    double pSum = 0;
    for(auto jet: jets) {
        if(!jet.passedJetSel()) continue;
	for(int i = 1; i <= 3; i++) {
	    sphere(i-1,i-1) += pow(jet[i], 2);
	    pSum += pow(jet[i], 2);
	    for(int j = i+1; j <= 3; j++) {
		sphere(i-1,j-1) += jet[i]*jet[j];
		sphere(j-1,i-1) += jet[i]*jet[j];
	    }
	}
    }
    sphere *= 1/pSum;
    TMatrixDSymEigen eigenGetter(sphere);
    TVectorD eigen = eigenGetter.GetEigenValues();
    return 3./2*(eigen(1) + eigen(2));
}

double JetCentrality(std::vector<GoodPart>& jets, double HT) {
    double eTot = 0;
    for(auto jet: jets) {
	if(!jet.passedJetSel()) continue;
	eTot += jet.E();
    }
    return HT/eTot;
}

std::pair<double,double> EventShape(std::vector<GoodPart>& jets, std::vector<GoodPart>& leps, double Met2, double MetPhi) {
    TMatrixDSym shape(3);
    double pSum = 0;
    for(auto jet: jets) {
	for(int i = 1; i <= 3; i++) {
	    shape(i-1,i-1) += pow(jet[i], 2);
	    pSum += pow(jet[i], 2);
	    for(int j = i+1; j <= 3; j++) {
		shape(i-1,j-1) += jet[i]*jet[j];
		shape(j-1,i-1) += jet[i]*jet[j];
	    }
	}
    }
    for(auto lep: leps) {
	for(int i = 1; i <= 3; i++) {
	    shape(i-1,i-1) += pow(lep[i], 2);
	    pSum += pow(lep[i], 2);
	    for(int j = i+1; j <= 3; j++) {
		shape(i-1,j-1) += lep[i]*lep[j];
		shape(j-1,i-1) += lep[i]*lep[j];
	    }
	}
    }
    pSum += Met2;
    shape(1, 1 ) += Met2*pow(cos(MetPhi), 2);
    shape(2, 2 ) += Met2*pow(sin(MetPhi), 2);
    shape(1, 2 ) += Met2*sin(MetPhi/2)/2;
    shape(2, 1 ) += Met2*sin(MetPhi/2)/2;
    shape *= 1/pSum;
    TMatrixDSymEigen eigenGetter(shape);
    TVectorD eigen = eigenGetter.GetEigenValues();
    return std::make_pair(3.*(eigen(1)*eigen(2) + eigen(1)*eigen(0) + eigen(0)*eigen(2)), 27*eigen(0)*eigen(1)*eigen(2));

}

#include "Math/GenVector/VectorUtil.h"

#endif /* HELPERFUNCTIONS_H */

