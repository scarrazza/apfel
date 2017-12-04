#pragma once
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <complex>

using namespace std;



namespace HELLx {


template <class T>
class vec2{

private:
  T a1, a2;

public:
  inline vec2(T b1=0., T b2=0.): a1(b1), a2(b2) {}
  inline vec2(const vec2 & a): a1(a.a1), a2(a.a2) {}

  inline T entry1() { return a1; };
  inline T entry2() { return a2; };

  //friend class sqmatrix;

  inline vec2 & operator = (const vec2 & a){
    if (a1 != a.a1 || a2 != a.a2) {
      a1 = a.a1; a2 = a.a2;
    }
    return *this;
  }
  inline vec2  operator - () const {
    vec2 A = *this;
    A.a1 = -a1; 
    A.a2 = -a2; 
    return A;
  }


  inline vec2& operator += (const vec2 &a) { a1 +=a.a1 ; a2 += a.a2; return *this; }
  inline vec2& operator -= (const vec2 &a) { a1 -=a.a1 ; a2 -= a.a2; return *this; }
  inline vec2& operator += (const T &a){ a1 += a; a2 += a; return *this;}
  inline vec2& operator -= (const T &a){ a1 -= a; a2 -= a; return *this;}
  inline vec2& operator *= (const T &a){ a1 *= a; a2 *= a; return *this;}
  inline vec2& operator /= (const T &a) { a1 /= a; a2 /= a; return *this; }

  //friend istream& operator >> (istream& is, vec2& a);
  //friend ostream& operator << (ostream& os,  const vec2& a);

  //friend vec2 operator + (const vec2 &a, const vec2 &b);
  //friend vec2 operator + (const vec2 &a, const T    &b);
  //friend vec2 operator + (const T    &b, const vec2 &a);

  //friend vec2 operator - (const vec2 &a, const vec2 &b);
  //friend vec2 operator - (const vec2 &a, const T    &b);
  //friend vec2 operator - (const T    &b, const vec2 &a);

  //friend vec2 operator * (const vec2 &a, const vec2 &b);
  //friend vec2 operator * (sqmatrix a, const vec2 &b);
  //friend vec2 operator * (const vec2 &a, T b);
  //friend vec2 operator * (T b, const vec2 &a);
  //friend vec2 operator / (const vec2 &a, const vec2 &b);
  //friend vec2 operator / (const vec2 &a, T b);
  //friend vec2 operator / (double b, const vec2 &a);
  //friend int operator == (const vec2 &a, const vec2 &b);
  //friend int operator != (const vec2 &a, const vec2 &b);

};




template <class T> vec2<T> operator + (const vec2<T> &a, const vec2<T> &b){
  vec2<T> c = a;  c += b;  return c;
}
template <class T> vec2<T> operator + (const vec2<T> &a, const T &b){
  vec2<T> c = a;  c += b;  return c;
}
template <class T> vec2<T> operator + (const T &b, const vec2<T> &a){
  vec2<T> c = a;  c += b;  return c;
}


template <class T> vec2<T> operator - (const vec2<T> &a, const vec2<T> &b){
  vec2<T> c = a;  c -= b;  return c;
}
template <class T> vec2<T> operator - (const vec2<T> &a, const T &b){
  vec2<T> c = a;  c -= b;  return c;
}
template <class T> vec2<T> operator - (const T &b, const vec2<T> &a){
  vec2<T> c = -a;  c += b;  return c;
}


template <class T> vec2<T> operator * (const vec2<T> &a, const T &b){
  vec2<T> c = a;  c *= b;  return c;
}
template <class T> vec2<T> operator * (const T &b, const vec2<T> &a){
  vec2<T> c = a;  c *= b;  return c;
}


template <class T> vec2<T> operator / (const vec2<T> &a, const T &b){
  vec2<T> c = a;  c /= b;  return c;
}

/*
int operator == (const vec2<T> &a, const vec2<T> &b) {
   return((a.r ==b.r) && (a.i == b.i) ? 1 : 0 );}

int operator != (const vec2<T> &a, const vec2<T> &b) {
   return((a.r == b.r) && (a.i == b.i) ? 0 : 1);}
*/

/*
template <class T> ostream & operator << (ostream & os, const vec2<T> & a){
  return os << "/ "  << setw(8) << setprecision(5) << a.entry1() << " \\" << endl
	    << "\\ " << setw(8) << setprecision(5) << a.entry2() << " /"  << endl;
}
*/
template <class T> ostream & operator << (ostream & os, vec2<T> & a){
  return os << "/ "  << setw(8) << setprecision(5) << a.entry1() << " \\" << endl
	    << "\\ " << setw(8) << setprecision(5) << a.entry2() << " /"  << endl;
}








template <class T>
class sqmatrix{

private:
  T a11, a12, a21, a22;

public:
  inline sqmatrix(T b11=0., T b12=0., T b21=0., T b22=0.): a11(b11), a12(b12), a21(b21), a22(b22) {}
  inline sqmatrix(const sqmatrix & a): a11(a.a11), a12(a.a12), a21(a.a21), a22(a.a22) {}

  inline T entry11() { return a11; }
  inline T entry12() { return a12; }
  inline T entry21() { return a21; }
  inline T entry22() { return a22; }

  inline T det() { return a11*a22 - a12*a21; }
  inline T trace() { return a11+a22; }

  inline sqmatrix & operator = (const sqmatrix & a){
    if (a11 != a.a11 || a12 != a.a12 || a21 != a.a21 || a22 != a.a22) {
      a11 = a.a11; a12 = a.a12; a21 = a.a21; a22 = a.a22;
    }
    return *this;
  }
  inline sqmatrix  operator - () const {
    sqmatrix A = *this;
    A.a11 = -a11; 
    A.a12 = -a12; 
    A.a21 = -a21; 
    A.a22 = -a22; 
    return A;
  }


  inline sqmatrix& operator += (const sqmatrix &a) {
    a11 += a.a11; a12 += a.a12; a21 += a.a21; a22 += a.a22; return *this;
  }
  inline sqmatrix& operator -= (const sqmatrix &a) {
    a11 -= a.a11; a12 -= a.a12; a21 -= a.a21; a22 -= a.a22; return *this;
  }
  inline sqmatrix& operator *= (const sqmatrix &b) {
    sqmatrix a = *this;
    a11 = a.a11*b.a11 + a.a12*b.a21;
    a12 = a.a11*b.a12 + a.a12*b.a22;
    a21 = a.a21*b.a11 + a.a22*b.a21;
    a22 = a.a21*b.a12 + a.a22*b.a22;
    return *this;
  }
  inline sqmatrix& operator += (const T &a) { a11 += a; a22 += a; return *this; }
  inline sqmatrix& operator -= (const T &a) { a11 -= a; a22 -= a; return *this; }
  inline sqmatrix& operator *= (const T &a) { a11 *= a; a12 *= a; a21 *= a; a22 *= a; return *this; }
  inline sqmatrix& operator /= (const T &a) { a11 /= a; a12 /= a; a21 /= a; a22 /= a; return *this; }

  inline sqmatrix transpose() {
    sqmatrix<T> b(a11, a21, a12, a22);
    return b;
  }

  inline sqmatrix inverse() {
    T d = det(*this);
    if(d==0) {
      cout << "Matrix not invertible!" << endl;
      exit(44);
    }
    sqmatrix<T> b(a22, -a12, -a21, a11);
    return b/d;
  }


  //friend istream& operator >> (istream& is, sqmatrix& a);
  //friend ostream& operator << (ostream& os,  const sqmatrix& a);

  //friend sqmatrix operator + (const sqmatrix &a, const sqmatrix &b);
  //friend sqmatrix operator + (const sqmatrix &a, T b);
  //friend sqmatrix operator + (T b, const sqmatrix &a);
  //
  //friend sqmatrix operator - (const sqmatrix &a, const sqmatrix &b);
  //friend sqmatrix operator - (const sqmatrix &a, T b);
  //friend sqmatrix operator - (T b, const sqmatrix &a);
  //friend sqmatrix operator * (const sqmatrix &a, const sqmatrix &b);
  //friend vec2<T>  operator * (const sqmatrix &a, vec2<T> b);
  //friend sqmatrix operator * (const sqmatrix &a, T b);
  //friend sqmatrix operator * (T b, const sqmatrix &a);
  //friend sqmatrix operator / (const sqmatrix &a, T b);
  //friend int operator == (const sqmatrix &a, const sqmatrix &b);
  //friend int operator != (const sqmatrix &a, const sqmatrix &b);


};






template <class T> inline sqmatrix<T> operator + (const sqmatrix<T> &a, const sqmatrix<T> &b) {
  sqmatrix<T> c = a;  c += b;  return c;
}
template <class T> inline sqmatrix<T> operator + (const sqmatrix<T> &a, const T &b) {
  sqmatrix<T> c = a;  c += b;  return c;
}
template <class T> inline sqmatrix<T> operator + (const T& b, const sqmatrix<T> &a){
  sqmatrix<T> c = a;  c += b;  return c;
}


template <class T> inline sqmatrix<T> operator - (const sqmatrix<T> &a, const sqmatrix<T> &b) {
  sqmatrix<T> c = a;  c -= b;  return c;
}
template <class T> inline sqmatrix<T> operator - (const sqmatrix<T> &a, const T &b) {
  sqmatrix<T> c = a;  c -= b;  return c;
}
template <class T> inline sqmatrix<T> operator - (const T& b, const sqmatrix<T> &a){
  sqmatrix<T> c = -a;  c += b;  return c;
}


template <class T> inline sqmatrix<T> operator * (const sqmatrix<T> &a, const sqmatrix<T> &b) {
  sqmatrix<T> c = a;  c *= b;  return c;
}
template <class T> inline sqmatrix<T> operator * (const sqmatrix<T> &a, const T &b) {
  sqmatrix<T> c = a;  c *= b;  return c;
}
template <class T> inline sqmatrix<T> operator * (const T& b, const sqmatrix<T> &a){
  sqmatrix<T> c = a;  c *= b;  return c;
}
template <class T, class U> inline sqmatrix<T> operator*(const sqmatrix<T> &a, const U &b) {
  sqmatrix<T> c = a;  c *= b;  return c;
}
template <class T, class U> inline sqmatrix<T> operator*(const U &b, const sqmatrix<T> &a) {
  sqmatrix<T> c = a;  c *= b;  return c;
}



template <class T> inline sqmatrix<T> operator / (const sqmatrix<T> &a, const T &b) {
  sqmatrix<T> c = a;  c /= b;  return c;
}
template <class T, class U> inline sqmatrix<T> operator / (const sqmatrix<T> &a, const U &b) {
  sqmatrix<T> c = a;  c /= b;  return c;
}





// template <class T> ostream & operator << (ostream & os, const sqmatrix<T> & a){
//   return os << "/ "  << setw(8) << setprecision(5) << a.entry11() << "  " << setw(8) << setprecision(5) << a.entry12() << " \\" << endl
// 	    << "\\ " << setw(8) << setprecision(5) << a.entry21() << "  " << setw(8) << setprecision(5) << a.entry22() << " /"  << endl;
// }

template <class T> ostream & operator << (ostream & os, sqmatrix<T> & a){
  return os << "/ "  << setw(8) << setprecision(5) << a.entry11() << "  " << setw(8) << setprecision(5) << a.entry12() << " \\" << endl
	    << "\\ " << setw(8) << setprecision(5) << a.entry21() << "  " << setw(8) << setprecision(5) << a.entry22() << " /"  << endl;
}




template <class T>
vec2<T> operator * (sqmatrix<T> &a, vec2<T> &b) {
  return vec2<T>(a.entry11()*b.entry1() + a.entry12()*b.entry2(),
		 a.entry21()*b.entry1() + a.entry22()*b.entry2() );
}


};


