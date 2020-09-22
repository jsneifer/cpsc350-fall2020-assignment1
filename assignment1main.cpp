#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <stdlib.h>

using namespace std;

void Generate(double m_implementMean, double m_implementVariance, double aStat, double cStat, double tStat, double gStat) {
  double y;
  double z;
  double m_standardGaussian;
  double m_randomD;
  int m_rand2;
  aStat *= 100; cStat *= 100; tStat *= 100; gStat *= 100;
  ofstream outputfile3;
  outputfile3.open("Sneifer.out", std::ios_base::app);
  outputfile3 << "\n1000 new lines: " << endl;
  if(outputfile3.is_open()) {
    for(int x = 0; x < 1000; ++x) {
      y = (double)rand() / RAND_MAX;
      z = (double)rand() / RAND_MAX;
      m_standardGaussian = sqrt(-2 * log(y)) * (cos(2 * M_PI * z));
      m_randomD = ((sqrt(m_implementVariance)) * m_standardGaussian) + m_implementMean;
        for(int q = 0; q < m_randomD; ++q) {
          m_rand2 = rand() % 100 + 1;
          if(m_rand2 > 1 && m_rand2 <= aStat) {
            outputfile3 << 'A';
          } else if(m_rand2 > aStat && m_rand2 <= (aStat+cStat)) {
            outputfile3 << 'C';
          } else if(m_rand2 > cStat && m_rand2 <= (aStat+cStat+tStat)) {
            outputfile3 << 'T';
          } else if(m_rand2 > tStat && m_rand2 <= (aStat+cStat+tStat+gStat)) {
            outputfile3 << 'G';
          }
        }
        outputfile3 << endl;
    }
  }
  outputfile3.close();
}

// Calulate return 1 if it is a valid text file
// or 0 if it is not a valid text file.
int Calculate(string textfile) {
  string tp;
  double m_mean = 0;
  int m_numLines = 0;
  int m_totalSize = 0;
  double m_variance = 0;
  double m_sdeviation;
  double a, c, t, g, aa, ac, at, ag, ca, cc, ct, cg, ta, tc, tt, tg, ga, gc, gt, gg=0;

  ifstream newFile;
  newFile.open(textfile);
  if(newFile.is_open()) {

    // Loops through the whole text file
    while(getline(newFile, tp)) {
      transform(tp.begin(), tp.end(), tp.begin(), ::toupper);
      ++m_numLines;
      m_mean += tp.size()-1;
      m_totalSize += tp.size()-1;

      // Loops through each character
      for(int i = 0; i < tp.size(); ++i) {
        if(tp[i] == 'A') {
          ++a;
          if(tp[i+1] == 'A') {
            ++aa;
            ++a;
          } else if(tp[i+1] == 'C') {
            ++ac;
            ++c;
          } else if(tp[i+1] == 'T') {
            ++at;
            ++t;
          }  else if(tp[i+1] == 'G') {
            ++ag;
            ++g;
          }
        } else if(tp[i] == 'C') {
          ++c;
          if(tp[i+1] == 'A') {
            ++ca;
            ++a;
          } else if(tp[i+1] == 'C') {
            ++cc;
            ++c;
          } else if(tp[i+1] == 'T') {
            ++ct;
            ++t;
          }  else if(tp[i+1] == 'G') {
            ++cg;
            ++g;
          }
        } else if(tp[i] == 'T') {
          ++t;
          if(tp[i+1] == 'A') {
            ++ta;
            ++a;
          } else if(tp[i+1] == 'C') {
            ++tc;
            ++c;
          } else if(tp[i+1] == 'T') {
            ++tt;
            ++t;
          }  else if(tp[i+1] == 'G') {
            ++tg;
            ++g;
          }
        } else if(tp[i] == 'G') {
          ++g;
          if(tp[i+1] == 'A') {
            ++ga;
            ++a;
          } else if(tp[i+1] == 'C') {
            ++gc;
            ++c;
          } else if(tp[i+1] == 'T') {
            ++gt;
            ++t;
          }  else if(tp[i+1] == 'G') {
            ++gg;
            ++g;
          }
        }
        ++i;
      }
    }
    newFile.close();

    ofstream outputfile2;
    outputfile2.open("Sneifer.out", std::ios_base::app);
    if(outputfile2.is_open()) {
      outputfile2 << "\nThe sum of the length of the DNA strings is: " << m_totalSize << endl;
    }

    // Probability of each NucleoTide
    a /= m_totalSize; c /= m_totalSize; t /= m_totalSize; g /= m_totalSize;
    // Probability of each bigragm
    m_totalSize /= 2;
    aa /= m_totalSize; ac /= m_totalSize; at /= m_totalSize; ag /= m_totalSize;
    ca /= m_totalSize; cc /= m_totalSize; ct /= m_totalSize; cg /= m_totalSize;
    ta /= m_totalSize; tc /= m_totalSize; tt /= m_totalSize; tg /= m_totalSize;
    ga /= m_totalSize; gc /= m_totalSize; gt /= m_totalSize; gg /= m_totalSize;

    // Mean calculation
    m_mean /= m_numLines;

    // Variance calculation
    ifstream sameFile(textfile);
    if(sameFile.is_open()) {
      while(getline(sameFile, tp)) {
        m_variance += pow(((tp.size()-1) - m_mean), 2);
      }
    }
    sameFile.close();
    m_variance /= m_numLines;
    // Standard Deviation Calculation
    m_sdeviation = sqrt(m_variance);

    if(outputfile2.is_open()) {
      outputfile2 << "Mean of lengths: " << m_mean<<endl;
      outputfile2 << "Variance of lengths: " << m_variance<<endl;
      outputfile2 << "Standard Deviation of lengths: " << m_sdeviation<<endl;
      outputfile2 << "\nRelative probability of nucleotides:"<<endl;
      outputfile2 << "A: "<<a<<"\nC: "<<c<<"\nT: "<<t<<"\nG: "<<g<<endl;
      outputfile2 << "\nRelative probability of nucleotide bigram: "<<endl;
      outputfile2 << "AA: "<<aa<<"\nAC: "<<ac<<"\nAT: "<<at<<"\nAG: "<<ag<<endl;
      outputfile2 << "CA: "<<ca<<"\nCC: "<<cc<<"\nCT: "<<ct<<"\nCG: "<<cg<<endl;
      outputfile2 << "TA: "<<ta<<"\nTC: "<<tc<<"\nTT: "<<tt<<"\nTG: "<<tg<<endl;
      outputfile2 << "GA: "<<ga<<"\nGC: "<<gc<<"\nGT: "<<gt<<"\nGG: "<<gg<<endl;

      outputfile2.close();
    }
  } else {
    cout << "Not a valid text file";
    return 0;
  }
  Generate(m_mean, m_variance, a, c, t, g);
  return 1;
}





int main(int argc, char** argv) {

  ofstream outputfile("Sneifer.out");
  if(outputfile.is_open()) {
    outputfile << "Joseph Sneifer\n2351513\nCPSC-350\n";
    outputfile.close();
  } else {
    cout << "Cannot open file";
  }


  Calculate(argv[1]);

  string rn;
  string yn;
  bool cont = true;
  while(cont == true) {
    cout << "\nWould you like to run another text file? (y/n)" << endl;
    getline(cin, yn);
    transform(yn.begin(), yn.end(), yn.begin(), ::tolower);
    if(yn == "y" || yn == "yes") {
      cout << "Enter a new text file: ";
      getline(cin, rn);
      Calculate(rn);
    } else if(yn == "n" || yn == "no") {
      cont == false;
      break;
    } else {
      cout << "That is not a valid response. " << endl;
    }
  }
  return 0;
}
