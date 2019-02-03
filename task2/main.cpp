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
using namespace std;


int main(int argc, const char * argv[]) {
    ifstream in;
    in.open("/Users/styskin/bio2019/task2/task2/0.txt");
    ofstream out;
    out.open("/Users/styskin/bio2019/task2/task2/0.out");
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
            cnk = cnk * (k - i + 1) / i;
            double p1 = pow(1.0 * (n) / L, i) * pow(1.0 * (L - n) / L, k - i);
            pp += cnk * p1;
            
            // two states and choosen majority
            // (i+1)/2 of i
            double px = 0;
            if (i % 2 == 0) {
                px += 0.5 * pow(p, i / 2) * pow(1 - p, i - i/2); // equal possibilities
            }
            // C N K
            for (int j = i/2 + 1; j <= i; ++j) {
                px += pow(p, j) * pow(1 - p, i - j); // all errors are errors
            }
            res += cnk * p1 * px;
        }
        cout << res * L << " prob = "  << pp  << endl;
    }
    in.close();
    out.close();
    
    return 0;
}
