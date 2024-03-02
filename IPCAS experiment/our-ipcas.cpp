extern "C" {
#include "miracl.h"
}
#include <cstdio>
#include <ctime>
#include "big.h"
#include "ecn.h"
#include <openssl/sha.h>
#include <sstream>
#include <vector>

Big randbits(int n) { Big z; bigbits(n, z.fn); return z; }
Miracl precision(196, 16);
miracl* mip = &precision;
ECn g, Pub;
Big q, alpha;
#define HASH_LEN 32


struct Sig
{
    ECn Y;
    Big w;
};

struct PK {
    ECn D;
    ECn R;
};

struct PSEU
{
    Big pid;
    long T;
};

void setup();
Big Hash(stringstream& st);

class Clas {
private:
public:
    virtual void reg() = 0;
    virtual Sig& sign(string& m) = 0;
    virtual PK& getPK() = 0;
    virtual ECn& getPub() = 0;
};

//! NIST p192 bits ECC curve prime
char* ecp = (char*)"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF";

//! NIST p192 bits ECC curve parameter b
char* ecb = (char*)"64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1";

//! NIST p192 bits ECC curve parameter q
char* ecq = (char*)"FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831";

//! NIST p192 bits ECC curve point of prime order (x,y)
char* ecx = (char*)"188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012";
char* ecy = (char*)"07192B95FFC8DA78631011ED6B24CDD573F977A11E794811";

class CLAS : public Clas {
public:
    virtual void reg() {
        Big x = rand(q),r = rand(q),d;
        X = x * g;
        pk.R = r * g;
        pseu = PSEU{ randbits(256), clock() };

        st << pseu.pid << pseu.T << pk.R << Pub;
        Big h2 = Hash(st);
        pk.D = h2 * X;
        pk.D += pk.R;

        st << pseu.pid << pseu.T << pk.D << pk.R << Pub;
        Big h1 = Hash(st);
        d = r + alpha * h1;

        sk = d + h2 * x;

    }

    PK& getPK() {
        return pk;
    }

    PSEU& getPSEU() {
        return pseu;
    }

    long getTimestamp() {
        return t;
    }

    ECn& getPub() {
        return Pub;
    }

    virtual Sig& sign(string& m) {
        Big y = rand(q);//random value
        ECn Y = y * g;
        t = clock();//timestamp

        st << m << pseu.pid << pseu.T << Y << pk.D << pk.R << t;
        Big h3 = Hash(st);

        sig.w = h3*sk+y;
        sig.Y = Y;
        return sig;
    }

private:
    PK pk;
    Big sk;
    Sig sig;
    ECn X;
    stringstream st;
    long t = 1;
    PSEU pseu;
};

void setup() {
    // Elliptic curve parameter reading
    Big a, b, p, px, py;
    int bits;
    a = -3;
    mip->IOBASE = 16;
    b = ecb;
    p = ecp;
    q = ecq;
    px = ecx;
    py = ecy;
    ecurve(a, b, p, MR_BEST);
    g = ECn(px, py);//generator
    alpha = rand(q);
    Pub = alpha * g;
}

Big Hash(stringstream& st) {
    size_t size = st.tellp();
    char* buff = new char[size];
    st.read(buff, size);
    unsigned char value[HASH_LEN];
    SHA256((unsigned char*)buff, size, value);
    st.str("");
    delete[] buff;
    return from_binary(HASH_LEN, (char*)value);
}

bool verify(Sig& sig, PSEU& pseu, PK& pk, ECn& Pub, string& m, long timestp) {

    stringstream st;
    ECn left, right;
    st << pseu.pid << pseu.T << pk.D << pk.R << Pub;
    Big h1 = Hash(st);
    st << m << pseu.pid << pseu.T << sig.Y << pk.D << pk.R << timestp;
    Big h3 = Hash(st);
    left = sig.w * g;
    left -= sig.Y;
    right = h1 * Pub;
    right += pk.D;
    right *= h3;
    if (left == right) {
        return true;
    }
    return false;
}

bool aggVerify(int n, string& msg, Big& aggSig, vector<PSEU>& vecPSEU, vector<PK>& vecPK, vector<ECn>& vecY, vector<long>& vecT) {
    ECn right;
    Big r1;
    stringstream st;
    ECn left = aggSig * g;
    for (int i = 0; i < n; i++) {
        st << vecPSEU[i].pid << vecPSEU[i].T << vecPK[i].D << vecPK[i].R << Pub;
        Big h1 = Hash(st);
        st << msg << vecPSEU[i].pid << vecPSEU[i].T << vecY[i] << vecPK[i].D << vecPK[i].R << vecT[i];
        Big h3 = Hash(st);
        
        left -= vecY[i];

        right += h3 * vecPK[i].D;
        r1 += h1 * h3;
    }
    right += r1 * Pub;
    if (left == right) {
        return true;
    }
    return false;
}

void singleTest(CLAS& clas) {
    double start;
    double diff;
    Sig sig;
    ECn left, right;
    string msg("This is a test.");
    cout << "First, we generate public key and scret key!" << endl;

    clas.reg();

    cout << "\nThen we compute the signature." << endl;
    start = clock();
    sig = clas.sign(msg);
    diff = ((double)clock() - start) / CLOCKS_PER_SEC * 1000.0;
    cout << "Sig: {" << sig.Y << ", " << sig.w << "}" << endl;
    printf("[*] Sign Time: %.6fms\n", diff);

    cout << "\nNow, we start to verify the sig." << endl;
    start = clock();
    if (verify(sig, clas.getPSEU(), clas.getPK(), clas.getPub(), msg, clas.getTimestamp())) {
        diff = ((double)clock() - start) / CLOCKS_PER_SEC * 1000.0;
        printf("[*] ACCEPT! Verification Time: %.6fms\n", diff);
    }
}

void avgTest(CLAS& clas, int n) {
    double s_start, v_start, s_total = 0, v_total = 0;
    Sig sig;
    ECn left, right;
    string msg("A traffic accident occurred 100 meters ahead.");
    for (int i = 0; i < n; i++) {
        clas.reg();
        s_start = clock();
        sig = clas.sign(msg);
        s_total += clock() - s_start;
        v_start = clock();
        if (verify(sig, clas.getPSEU(), clas.getPK(), clas.getPub(), msg, clas.getTimestamp())) {
            v_total += clock() - v_start;
        }
        else {
            cout << "[x] verification reject!" << endl;
            exit(-1);
        }
    }
    printf("[*] Average Signing Time: %.6fms\n", s_total / n / CLOCKS_PER_SEC * 1000.0);
    printf("[*] Average Individual Verification Time: %.6fms\n", v_total / n / CLOCKS_PER_SEC * 1000.0);
}

void aggTest(CLAS& ipcas, int n) {
    string msg("This is a test.");
    vector<ECn> vecY;
    vector<PSEU> vecPSEU;
    vector<PK> vecPK;
    vector<long> vecT;
    stringstream st;
    Big aggSig(0);
    for (int i = 0; i < n; i++) {
        ipcas.reg();
        Sig sig = ipcas.sign(msg);
        st << msg << ipcas.getPSEU().pid << ipcas.getPSEU().T << ipcas.getPK().D << ipcas.getPK().R << ipcas.getTimestamp();
        Big h4 = Hash(st);
        aggSig += sig.w;
        vecY.push_back(sig.Y);
        vecPSEU.push_back(ipcas.getPSEU());
        vecPK.push_back(ipcas.getPK());
        vecT.push_back(ipcas.getTimestamp());
    }
    double start = clock();
    if (aggVerify(n, msg, aggSig, vecPSEU, vecPK, vecY, vecT)) {
        double end = clock();
        printf("[*] %d Aggregate Verification Time: %.6fms\n", n, (end - start) / CLOCKS_PER_SEC * 1000.0);
    }
}

int main() {
    irand(2022l); // Set random seeds
    setup();
    cout << endl << "------------IPCAS-------------" << endl;
    CLAS ipcas;
    //singleTest(ipcas);
    avgTest(ipcas, 1000);//average signle signature&verfication excution time
    aggTest(ipcas, 20);//average aggregate signature verfication excution time
    aggTest(ipcas, 40);
    aggTest(ipcas, 60);
    aggTest(ipcas, 80);
    aggTest(ipcas, 100);
    return 0;
}