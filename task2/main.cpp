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
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace std;


struct TNuc {
    int first;
    int second;
    int third;
    int fourth;
    
    TNuc(int a, int b, int c, int d) {
        first = a;
        second = b;
        third = c;
        fourth = d;
    }
    
    TNuc(const TNuc& nuc) {
        first = nuc.first;
        second = nuc.second;
        third = nuc.third;
        fourth = nuc.fourth;
    }
    
    bool operator==(const TNuc &other) const {
        return (first == other.first
              && second == other.second
              && third == other.third);
    }
};

template <>
struct hash<TNuc>
{
    size_t operator()( const TNuc& k ) const
    {
        // Compute individual hash values for first, second and third
        // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
        res = res * 31 + hash<int>()( k.first );
        res = res * 31 + hash<int>()( k.second );
        res = res * 31 + hash<int>()( k.third );
        res = res * 31 + hash<int>()( k.fourth );
        return res;
    }
};

typedef unordered_map<TNuc, double> TMap;


int main(int argc, const char * argv[]) {
    ifstream in;
    in.open("/Users/styskin/bio2019/task2/task2/2.txt");
    ofstream out;
    out.open("/Users/styskin/bio2019/task2/task2/2.out");
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
            if (0) { // easy problem
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
            } else {  // hard problem
                vector<TMap> M(2);
                M[0][TNuc(1, 0, 0, 0)] = (1-p);
                M[0][TNuc(0, 1, 0, 0)] = p / 3;
                M[0][TNuc(0, 0, 1, 0)] = p / 3;
                M[0][TNuc(0, 0, 0, 1)] = p / 3;

                for (int j = 0; j < i-1; ++j) {
                    M[(j+1) % 2].clear();
                    for (const auto& nuc : M[j%2]) {
                        TNuc k1(nuc.first);
                        TNuc k2(nuc.first);
                        TNuc k3(nuc.first);
                        TNuc k4(nuc.first);
                        ++k1.first;
                        ++k2.second;
                        ++k3.third;
                        ++k4.fourth;
                        M[(j+1) % 2][k1] += nuc.second * (1-p);
                        M[(j+1) % 2][k2] += nuc.second * p / 3;
                        M[(j+1) % 2][k3] += nuc.second * p / 3;
                        M[(j+1) % 2][k4] += nuc.second * p / 3;
                    }
                }
                for (const auto& nuc : M[(i-1) % 2]) {
                    int mx = max(nuc.first.first, max(nuc.first.second, max(nuc.first.third, nuc.first.fourth)));
                    int c = 0;
                    if (nuc.first.first == mx) {
                        ++c;
                    }
                    if (nuc.first.second == mx) {
                        ++c;
                    }
                    if (nuc.first.third == mx) {
                        ++c;
                    }
                    if (nuc.first.fourth == mx) {
                        ++c;
                    }
                    if (nuc.first.first == mx && c > 1) {
                        res += p1 * nuc.second * (c - 1) / c;
                    } else if (nuc.first.first < mx) {
                        res += p1 * nuc.second;
                    }
                }
                
            }
        }
        out << setprecision(16) << res * L << endl;
    }
    in.close();
    out.close();
    
    return 0;
}
