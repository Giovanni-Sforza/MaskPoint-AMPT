#ifndef CORRELATION_CALCULATOR_H
#define CORRELATION_CALCULATOR_H

#include <iostream>
#include "TComplex.h"
#include "TObject.h"

class CorrelationCalculator : public TObject 
{
    public:
    CorrelationCalculator(){}
    // static const int 20 = 20;
    // static const int 20 = 20; 

    TComplex Qvector[20][20];
    TComplex Qvector0M[20][20];
    TComplex Qvector0P[20][20];
    TComplex Qvector2M[20][20];
    TComplex Qvector2P[20][20];
    TComplex Qvector4M[20][20];
    TComplex Qvector4P[20][20];
    TComplex Qvector6M[20][20];
    TComplex Qvector6P[20][20];
    TComplex Qvector8M[20][20];
    TComplex Qvector8P[20][20];
    TComplex Qvector10M[20][20];
    TComplex Qvector10P[20][20];
    TComplex Qvector12M[20][20];
    TComplex Qvector12P[20][20];
    TComplex Qvector14M[20][20];
    TComplex Qvector14P[20][20];
    TComplex Qvector16M[20][20];
    TComplex Qvector16P[20][20];
    TComplex Qvector18M[20][20];
    TComplex Qvector18P[20][20];
    TComplex QvectorSubLeft[20][20];
    TComplex QvectorSubMiddle[20][20];
    TComplex QvectorSubRight[20][20];
    TComplex pvector[20][20];
    TComplex pvectorM[20][20];
    TComplex pvectorP[20][20];
    TComplex qvector[20][20];
    TComplex pvector0M[20][20];
    TComplex pvector0P[20][20];
    TComplex pvector4M[20][20];
    TComplex pvector4P[20][20];
    TComplex pvector8M[20][20];
    TComplex pvector8P[20][20];

    TComplex Two(int n1, int n2);
    TComplex Two_3SubLM(int n1, int n2);
    TComplex Two_3SubRM(int n1, int n2);
    TComplex Two_3SubLR(int n1, int n2);

    TComplex Three(int n1, int n2, int n3);
    TComplex Three_3Sub(int n1, int n2, int n3);

    TComplex Four(int n1, int n2, int n3, int n4);
    TComplex Four_3SubMMLR(int n1, int n2, int n3, int n4);
    TComplex Four_3SubLLMR(int n1, int n2, int n3, int n4);
    TComplex Four_3SubRRML(int n1, int n2, int n3, int n4);
    TComplex Five(int n1, int n2, int n3, int n4, int n5);
    TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
    TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
    TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);

    TComplex Q(int n, int p);
    TComplex QsubLeft(int n, int p);
    TComplex QsubMiddle(int n, int p);
    TComplex QsubRight(int n, int p);
    TComplex p(int n, int p);
    TComplex q(int n, int p);

    void ResetQ(const int nMaxHarm, const int nMaxPow);
    void FillQVector(TComplex _Qvector[20][20], double _Qcos[20][20], double _Qsin[20][20]);

    protected:
    private:
    // ClassDef(CorrelationCalculator, 1);   //
};
#endif
