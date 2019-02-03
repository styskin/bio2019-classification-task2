//
//  main.cpp
//  task2
//
//  Created by Andrey Styskin on 02.02.2019.
//  Copyright © 2019 Andrey Styskin. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;


int main(int argc, const char * argv[]) {
    ifstream in;
    in.open("/Users/styskin/bio2019/task2/task2/1.txt");
    ofstream out;
    out.open("/Users/styskin/bio2019/task2/task2/1.out");
    int T;
    in >> T;
    for (int t = 0; t < T; ++t) {
        int L, n, k;
        double p;
        in >> L >> n >> p >> k;
        // L - length
        // n - read length
        // p
        // k - number of reads
        double res = 0;
        double cnk = 1;
        double pp = 0;
        double tp = pow(1.0 * (L - n) / L, k);
        pp += tp;
        res += tp * 0.75;
        for (int i = 1; i <= k; ++i) {
            // 0.214491825
            // 0,034416825
            cnk = cnk * (k - i + 1) / i;
            double p1 = cnk * pow(1.0 * (n) / L, i) * pow(1.0 * (L - n) / L, k - i);
            pp += p1;

            double px = 0;
            vector< vector<double> > M(i, vector< double >(i + 1));

//            vector< vector<double> > M(i, vector< double >(i));  // i x i matrix
            M[0][0] = (1-p);
            M[0][1] = p;
            for (int j = 0; j < i-1; ++j) {
                for (int k = 0; k < i; ++k) {
                    M[j+1][k]     += M[j][k] * (1 - p);
                    M[j+1][k + 1] += M[j][k] * p;
                }
            }
//            double CHECK = 0;
//            for (int k = 0; k <= i; ++k) {
//                CHECK += M[i-1][k];
//            }
//            cout << "Check " << CHECK << endl;
            if (i % 2 == 0) {
                res += 0.5 * p1 * M[i-1][i / 2];
            }
            for (int k = i / 2 + 1; k <= i; ++k) {
                res += p1 * M[i-1][k];
            }
            
    
//            double cnk2 = 1;
//            for (int j = 1; j < i / 2; ++j) {
//                cnk2 = cnk2 * (i - j + 1) / j;
//            }
//            if (i % 2 == 0) { // equal probabilities
//                cnk2 = cnk2 * (i - (i/2) + 1) / (i/2);
//                px += 0.5 * pow(p, i / 2) * pow(1 - p, i - i/2) * cnk2;
//            }
//            // C N K
//            for (int j = i/2 + 1; j <= i; ++j) {
//                cnk2 = cnk2 * (i - j + 1) / j;
//                px += pow(p, j) * pow(1 - p, i - j) * cnk2; // all errors are errors
//            }
//            res += p1 * px;
        }
        out << res * L << endl;
    }
    in.close();
    out.close();
    
    return 0;
}
