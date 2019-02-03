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
        double pp = 0;
        double tp = pow(1.0 * (L - n) / L, k);
        pp += tp;
        res += tp * 0.75;
        double p1 = pow(1.0 * (L - n) / L, k);
        for (int i = 1; i <= k; ++i) {
            p1 *= (1.0 * (k - i + 1) / i) * (1.0 * n / (L - n));
            pp += p1;
            vector< vector<double> > M(2, vector< double >(i + 1));
            M[0][0] = (1-p);
            M[0][1] = p;
            for (int j = 0; j < i-1; ++j) {
                M[(j+1) % 2][0] = 0;
                for (int k = 0; k < i; ++k) {
                    M[(j+1) % 2][k]     += M[j%2][k] * (1 - p);
                    M[(j+1) % 2][k + 1] = M[j%2][k] * p;
                }
            }
            if (i % 2 == 0) {
                res += 0.5 * p1 * M[(i-1) % 2][i / 2];
            }
            for (int k = i / 2 + 1; k <= i; ++k) {
                res += p1 * M[(i-1) % 2][k];
            }
        }
//        cout << "Check " << pp << endl;
        out << res * L << endl;
    }
    in.close();
    out.close();
    
    return 0;
}
