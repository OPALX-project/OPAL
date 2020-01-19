// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

/***************************************************************************

 TEST 1: testing accumulation: Field += Field
 TEST 2: testing substraction: Field -= Field
 TEST 3: testing multipplication: Field *= Field
 TEST 4: testing division: Field *= Field

 TEST 5: testing scalar addition: Field += T
 TEST 6: testing scalar substraction: Field -= T
 TEST 7: testing scalar multipplication: Field *= T
 TEST 8: testing scalar division: Field /= T

 TEST 9: testing field expression and accumulation: Field += ExprElem
 TEST 10: testing field expression and substraction: Field -= ExprElem
 TEST 11: testing field expression and multipplication: Field *= ExprElem
 TEST 12: testing field expression and division: Field /= ExprElem

 TEST 13: testing accumulation with indexed field: IndexingField += IndexingField
 TEST 14: testing substraction with indexed field: IndexingField -= IndexingField
 TEST 15: testing multipplication with indexed field: IndexingField *= IndexingField
 TEST 16: testing division with indexed field: IndexingField /= IndexingField

 TEST 17: testing accumulation with indexed field: IndexingField += T
 TEST 18: testing substraction with indexed field: IndexingField -= T
 TEST 19: testing multipplication with indexed field: IndexingField *= T
 TEST 20: testing division with indexed field: IndexingField /= T

 TEST 21: testing accumulation index: IndexingField += Index
 TEST 22: testing substraction index: IndexingField -= Index
 TEST 23: testing multipplication index: IndexingField *= Index
 TEST 24: testing division index: IndexingField /= Index

 TEST 25: testing accumulation with IndexingFields on rhs: IndexingField += ExprElem
 TEST 26: testing substraction with IndexingFields on rhs: IndexingField -= ExprElem
 TEST 27: testing multipplication with IndexingFields on rhs: IndexingField *= ExprElem
 TEST 28: testing divisiond with IndexingFields on rhs: IndexingField /= ExprElem

 // testing accumulation: IndexingField += ExprElem
 // (with Fields on rhs

 //  // test 29
 //  A = 0.0;
 //  B = 1.0;
 //  C = 2.0;
 //  A[I][J] += B + C;
 //  testmsg << A << endl;
 usage: chsr-2 <gridSize>

 note: gridSize is used for all 6 dimensions in this testprogram

***************************************************************************/

#include "Ippl.h"

// dimension of f
const unsigned Dim = 6;

// presision
typedef double T;

// typedef our particle layout type
typedef ParticleSpatialLayout<T,Dim> playout_t;

int n;

const T pi = acos(-1.0);
const T fMax = 1.0;       // maximum value for f

void notify( bool passedSingleTest, int *testNum);


void notify( bool passedSingleTest, int *testNum, Timer *t) {
  Inform msg;

  if (!passedSingleTest)
    msg << "test "<<  *testNum << " failed" << endl;
  else {
    msg << "test "<<  *testNum << " passed";
    msg << " CpuTime= " << t->clock_time() << endl;
  }
  (*testNum)++;
}



int main(int argc, char *argv[]){
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);

  n = atoi(argv[1]);

  msg << "Dim=" << Dim << " n= "<< n << " M= " << pow(Dim,n) << endl;

  // create layout objects for max DIM==6
  const Index I(n), J(n), K(n);          // spacial
  const Index L(n), N(n), M(n);          // momenta

  // Initialize domain objects

  NDIndex<Dim> domain, domain1;

  domain[0] = n;
  domain[1] = n;
  domain[2] = n;
  domain[3] = n;

  domain1[0] = n+1;
  domain1[1] = n+1;
  domain1[2] = n+1;
  domain1[3] = n+1;

  if (Dim==6) {
    domain[4] = n;
    domain[5] = n;
    domain1[4] = n+1;
    domain1[5] = n+1;
  }

  UniformCartesian<Dim> mymesh(domain1);
  FieldLayout<Dim> FL(domain);

  int testNum = 1;
  double eps = 1.0e-07;

  // Flags to keep track of pass/fail results:
  bool passed = true;           // True at end means passed all tests
  bool passedSingleTest = true; // For individual tests

  Timer timer1;

  timer1.clear();
  timer1.start();
  Field<T,Dim> A(mymesh,FL);
  Field<T,Dim> B(mymesh,FL);
  Field<T,Dim> C(mymesh,FL);
  timer1.stop();
  msg << "Static field objects created t= " << timer1.clock_time() << endl;

  Field<double,Dim>::iterator fIt;
  T a,b,c;


  // TEST 1: testing accumulation: Field += Field
  a = 0.0;   b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;   B = 1.0;
  A += B;
  timer1.stop();

  a += b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 2: testing substraction: Field -= Field
  a = 0.0;   b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;   B = 1.0;
  A -= B;
  timer1.stop();

  a -= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 3: testing multipplication: Field *= Field
  a = 1.0;   b = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;   B = 2.0;
  A *= B;
  timer1.stop();

  a *= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 4: testing multipplication: Field /= Field
  a = 1.0;   b = 5.0;
  timer1.clear();
  timer1.start();
  A = 1.0;   B = 5.0;
  A /= B;
  timer1.stop();

  a /= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 5: testing scalar addition: Field += T
  a = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  A += a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 6: testing scalar substraction: Field -= T
  a = 1.0;
  timer1.clear();
  timer1.start();
  A = 2.0;
  A -= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 7: testing scalar multipplication: Field *= T
  a = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  A *= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 8: testing scalar division: Field /= T
  a = 2.0;
  timer1.clear();
  timer1.start();
  A = 4.0;
  A /= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 9: testing field expression and accumulation: Field += ExprElem
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  C = 2.0;
  A += B + C;
  timer1.stop();

  a += b + c;

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 10: testing field expression and substraction: Field -= ExprElem
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  C = 2.0;
  A -= B + C;
  timer1.stop();

  a -= b + c;

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 11: testing field expression and multipplication: Field *= ExprElem
  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  C = 3.0;
  A *= B + C;
  timer1.stop();

  a *= b + c;

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 12: testing field expression and division: Field /= ExprElem

  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  C = 3.0;
  A /= B + C;
  timer1.stop();

  a /= b + c;

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 13: testing accumulation with indexed field: IndexingField += IndexingField
  a = 0.0;
  b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  A[I][J][K][L][M][N] += B[I][J][K][L][M][N];
  timer1.stop();
  a += b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 14: testing substraction with indexed field: IndexingField -= IndexingField
  a = 0.0;
  b = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  B = 1.0;
  A[I][J][K][L][M][N] -= B[I][J][K][L][M][N];
  timer1.stop();
  a -= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 15: testing multipplication with indexed field: IndexingField *= IndexingField
  a = 1.0;
  b = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  A[I][J][K][L][M][N] *= B[I][J][K][L][M][N];
  timer1.stop();
  a *= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 16: testing division with indexed field: IndexingField /= IndexingField
  a = 1.0;
  b = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  B = 2.0;
  A[I][J][K][L][M][N] /= B[I][J][K][L][M][N];
  timer1.stop();
  a /= b;
  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 17: testing accumulation with indexed field: IndexingField += T
  a = 1.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  A[I][J][K][L][M][N] += a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 18: testing substraction with indexed field: IndexingField -= T
  a = 1.0;
  timer1.clear();
  timer1.start();
  A = 2.0;
  A[I][J][K][L][M][N] -= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 19: testing multipplication with indexed field: IndexingField *= T
  a = 2.0;
  timer1.clear();
  timer1.start();
  A = 1.0;
  A[I][J][K][L][M][N] *= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 20: testing division with indexed field: IndexingField /= T
  a = 2.0;
  timer1.clear();
  timer1.start();
  A = 4.0;
  A[I][J][K][L][M][N] /= a;
  timer1.stop();

  passedSingleTest = true;

  for (fIt=A.begin() ; fIt!=A.end(); ++fIt){
    if (fabs(*fIt - a) > eps) {
      passedSingleTest = false; passed = false;
      fIt = A.end();
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 21: testing accumulation index: IndexingField += Index
  a = 50.0;
  timer1.clear();
  timer1.start();
  A = 0.0;
  A[I][J][0][0][0][0] += I;
  timer1.stop();
  T sumval = sum(A);
  passedSingleTest = true;

  if (fabs(sumval - a) > eps) {
    passedSingleTest = false; passed = false;
    fIt = A.end();
  }

  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 22: testing substraction index: IndexingField -= Index
  // TEST 23: testing multipplication index: IndexingField *= Index
  // TEST 24: testing division index: IndexingField /= Index

  testNum = 25;

  /*
    Note: have to use reduce in order to deduce the result
          of the "local" compares of the field values
  */

  NDIndex<Dim> idx = FL.getLocalNDIndex();
  NDIndex<Dim> elem;

  // TEST 25: testing accumulation with IndexingFields on rhs: IndexingField += ExprElem
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] += B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();

  a += b + c;

  passedSingleTest = true;

  for (int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
        passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 26: testing substraction with IndexingFields on rhs: IndexingField -= ExprElem
  a = 0.0;
  b = 1.0;
  c = 2.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] -= B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();

  a -= b + c;

  passedSingleTest = true;

  for (int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
        passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------


  // TEST 27: testing multipplication with IndexingFields on rhs: IndexingField *= ExprElem
  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] *= B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();

  a *= b + c;

  passedSingleTest = true;

  for (int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
        passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  // TEST 28: testing divisiond with IndexingFields on rhs: IndexingField /= ExprElem
  a = 1.0;
  b = 2.0;
  c = 3.0;
  timer1.clear();
  timer1.start();
  A = a;
  B = b;
  C = c;
  A[I][J][0][0][0][0] /= B[I][J][0][0][0][0] + C[I][J][0][0][0][0];
  timer1.stop();

  a /= b + c;

  passedSingleTest = true;

  for (int i=idx[0].min(); i < idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (int j=idx[1].min(); j < idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      if (fabs(A.localElement(elem) - a) > eps) {
        passedSingleTest = false; passed = false;
      }
    }
  }
  notify(passedSingleTest,&testNum,&timer1);
  // --------------------------------------------

  Ippl::Comm->barrier();
  msg << "done ...." << endl;
  return 0;
}
