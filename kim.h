
#define IDX(i,j) ( (i) + n * (j) )
void derivKim4_PQ(double *P, double *Q, int n)
{
    const double alpha = 0.5862704032801503;
    const double beta = 9.549533555017055e-2;

    const double a1 = 0.6431406736919156;
    const double a2 = 0.2586011023495066;
    const double a3 = 7.140953479797375e-3;

    const double y00 = 0.0;
    const double y10=8.360703307833438e-2;
    const double y20=3.250008295108466e-2;
    const double y01=5.912678614078549;
    const double y11= 0.0;
    const double y21=0.3998040493524358;
    const double y02=3.775623951744012;
    const double y12=2.058102869495757;
    const double y22=0.0;
    const double y03=0.0;
    const double y13=0.9704052014790193;
    const double y23=0.7719261277615860;
    const double y04=0.0;
    const double y14=0.0;
    const double y24=0.1626635931256900;

    const double b10=-0.3177447290722621;
    const double b20=-0.1219006056449124;
    const double b01=-3.456878182643609;
    const double b21=-0.6301651351188667;
    const double b02=5.839043358834730;
    const double b12=-2.807631929593225e-2;
    const double b03=1.015886726041007;
    const double b13=1.593461635747659;
    const double b23=0.6521195063966084;
    const double b04=-0.2246526470654333;
    const double b14=0.2533027046976367;
    const double b24=0.3938843551210350;
    const double b05=8.564940889936562e-2;
    const double b15=-3.619652460174756e-2;
    const double b25=1.904944407973912e-2;
    const double b06=-1.836710059356763e-2;
    const double b16=4.080281419108407e-3;
    const double b26=-1.027260523947668e-3;

    const double b00 = -(b01 + b02 + b03 + b04 + b05 + b06);
    const double b11 = -(b10 + b12 + b13 + b14 + b15 + b16);
    const double b22 = -(b20 + b21 + b23 + b24 + b25 + b26);

    const int nd = n * n;

    for (int i = 0; i < nd; i++) {
        P[i] = 0.0;
    }
    for (int i = 3; i < n-3; i++) {
        P[IDX(i,i-2)] = beta;
        P[IDX(i,i-1)] = alpha;
        P[IDX(i,i)] = 1.0;
        P[IDX(i,i+1)] = alpha;
        P[IDX(i,i+2)] = beta;
    }


    P[IDX(0,0)] = 1.0;
    P[IDX(0,1)] = y01;
    P[IDX(0,2)] = y02;

    P[IDX(1,0)] = y10;
    P[IDX(1,1)] = 1.0;
    P[IDX(1,2)] = y12;
    P[IDX(1,3)] = y13;

    P[IDX(2,0)] = y20;
    P[IDX(2,1)] = y21;
    P[IDX(2,2)] = 1.0;
    P[IDX(2,3)] = y23;
    P[IDX(2,4)] = y24;

    P[IDX(n-3,n-5)] = y24;
    P[IDX(n-3,n-4)] = y23;
    P[IDX(n-3,n-3)] = 1.0;
    P[IDX(n-3,n-2)] = y21;
    P[IDX(n-3,n-1)] = y20;

    P[IDX(n-2,n-4)] = y13;
    P[IDX(n-2,n-3)] = y12;
    P[IDX(n-2,n-2)] = 1.0;
    P[IDX(n-2,n-1)] = y10;

    P[IDX(n-1,n-3)] = y02;
    P[IDX(n-1,n-2)] = y01;
    P[IDX(n-1,n-1)] = 1.0;

    for (int j = 0; j < nd; j++) {
        Q[j] = 0.0;
    }
    for (int i = 3; i < n-3; i++) {
        Q[IDX(i,i-3)] = -a3;
        Q[IDX(i,i-2)] = -a2;
        Q[IDX(i,i-1)] = -a1;
        Q[IDX(i,i)] = 0.0;
        Q[IDX(i,i+1)] = a1;
        Q[IDX(i,i+2)] = a2;
        Q[IDX(i,i+3)] = a3;
    }

    Q[IDX(0,0)] = b00;
    Q[IDX(0,1)] = b01;
    Q[IDX(0,2)] = b02;
    Q[IDX(0,3)] = b03;
    Q[IDX(0,4)] = b04;
    Q[IDX(0,5)] = b05;
    Q[IDX(0,6)] = b06;

    Q[IDX(1,0)] = b10;
    Q[IDX(1,1)] = b11;
    Q[IDX(1,2)] = b12;
    Q[IDX(1,3)] = b13;
    Q[IDX(1,4)] = b14;
    Q[IDX(1,5)] = b15;
    Q[IDX(1,6)] = b16;

    Q[IDX(2,0)] = b20;
    Q[IDX(2,1)] = b21;
    Q[IDX(2,2)] = b22;
    Q[IDX(2,3)] = b23;
    Q[IDX(2,4)] = b24;
    Q[IDX(2,5)] = b25;
    Q[IDX(2,6)] = b26;

    Q[IDX(n-3,n-1)]   = -b20;
    Q[IDX(n-3,n-2)] = -b21;
    Q[IDX(n-3,n-3)] = -b22;
    Q[IDX(n-3,n-4)] = -b23;
    Q[IDX(n-3,n-5)] = -b24;
    Q[IDX(n-3,n-6)] = -b25;
    Q[IDX(n-3,n-7)] = -b26;

    Q[IDX(n-2,n-1)]   = -b10;
    Q[IDX(n-2,n-2)] = -b11;
    Q[IDX(n-2,n-3)] = -b12;
    Q[IDX(n-2,n-4)] = -b13;
    Q[IDX(n-2,n-5)] = -b14;
    Q[IDX(n-2,n-6)] = -b15;
    Q[IDX(n-2,n-7)] = -b16;

    Q[IDX(n-1,n-1)]   = -b00;
    Q[IDX(n-1,n-2)] = -b01;
    Q[IDX(n-1,n-3)] = -b02;
    Q[IDX(n-1,n-4)] = -b03;
    Q[IDX(n-1,n-5)] = -b04;
    Q[IDX(n-1,n-6)] = -b05;
    Q[IDX(n-1,n-7)] = -b06;

}
