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
        for (int i = 0; i < k; ++i) {
            if (i == 0) {
                res += pow(1.0 * (L - n) / L, k) * 0.75;
            } else {
                double p1 = pow(1.0 * (n) / L, i);
                // two states and choosen majority
                // (i+1)/2 of i
                double px = 0;
                if (i % 2 == 0) {
                    px += 0.5 * pow(p, i / 2); // equal possibilities
                }
                for (int j = (i + 3)/2; j <= i; ++j) {
                    px += pow(p, j); // all errors are errors
                }
                res += p1 * px;
            }
        }
        out << res * L << endl;
    }
    in.close();
    out.close();
    
    return 0;
}
