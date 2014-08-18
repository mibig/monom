#include <iostream>
#include <vector>

int regul(std::vector<int> v) {
  while (v.vector::back() == 0) {
    v.vector::pop_back();
  }
  return 0;
} 

class zzbar{

  public: 
    zzbar(int d, int n) {
      if (d > 1 && n <d) {
        deg = d;
        num = n;
      }
    }
    // set deg and num
    void set_dg(int d, int n) {
      deg = d;
      num = n;
    }
    // Function finding \alpha          
    int alpha() {
      int alp, halfnum;
      halfnum = (1 + num) / 2;
      alp = deg - halfnum;   
      return alp;
    }
    // Function finding \beta           
    int beta() {
      int bet;
      bet = (1 + num) / 2;
      return bet;
    }
    // Function finding prefix          
    int prefix() {
      int prefix;
      if (deg%2==0 && num==deg-1) {
                     prefix = 3;    // without prefix          
      } 
      else {
        if(num%2==1) prefix = 1;    // prefix = Re
        if(num%2==0) prefix = 2;    // prefix = Im
      } 
      return prefix;
    }
    // Operator < comparing two monoms zzbar
    bool operator<(zzbar ob2) {
      bool answ;
      if (deg < ob2.deg) {
        answ = true;
      }
      else if (deg == ob2.deg && num > ob2.num) {
        answ = true;
      }
      else answ = false;

      return answ;
    }
    // Operator == comparing two monoms zzbar
    bool operator==(zzbar ob2) {
      bool answ;
      if (deg == ob2.deg && num == ob2.num) {
        answ = true;
      }
      else answ = false;

      return answ;
    }
    // Function printing monom zzbar          
    void print_zzbar() {
      int alp, bet, pr;
      pr  = prefix();
      alp = alpha();
      bet = beta();
      if (pr == 1) {
        std::cout << "\\mathrm{Re}\\left(";
      }
      else if (pr == 2) {
        std::cout << "\\mathrm{Im}\\left(";
      }
      if (alp == 1) {
        std::cout << "z";
      }
      else {
        std::cout << "z^{" << alp << "}";
      }
      if (bet == 1) {
        std::cout << "\\bar{z}";
      }
      else {
        std::cout << "\\bar{z}^{" << bet << "}";
      }
      if(pr < 3) {
        std::cout << "\\right)";
      }
    }
  private:
    int deg;  // degree of monom z^{\alpha}\bar{z}^{\beta}
    int num;  // numer of monom among all monoms degree deg
};

class monom : public zzbar{ 
  public: 
    monom(int d, int n, std::vector<int> vec):zzbar(d, n) {
      degu = vec;
    }
/*    // Operator < comparing two monoms zzbar
    bool operator<(monom ob2) {
      bool answ;
      if (deg < ob2.deg) {
        answ = true;
      }
      else if (deg == ob2.deg && num > ob2.num) {
        answ = true;
      }
      else answ = false;

      return answ;
    }
    // Operator == comparing two monoms zzbar
    bool operator==(zzbar ob2) {
      bool answ;
      if (deg == ob2.deg && num == ob2.num) {
        
        if () {
          
        }
        else answ = false;
      }
      else answ = false;

      return answ;
    }   */
    // Function printing monom          
    void print_monom() {
      int i, n = degu.vector::size();
      print_zzbar();
      for (i=0; i<n; ++i) {
        if (degu[i] > 1) {
          std::cout << "u_{" << i+1 << "}^{" << degu[i] << "}";
        }
        else if (degu[i] == 1) {
          std::cout << "u_{" << i+1 << "}";
        }
      }
      std::cout << "\n";
    }
  private: 
    std::vector<int> degu;   
};

int main() 
{ 

  std::vector<int> vect;
  vect.vector::push_back(1);
  vect.vector::push_back(7);
  vect.vector::push_back(4);
  vect.vector::push_back(0);
  vect.vector::push_back(12);
  vect.vector::push_back(0);
  vect.vector::push_back(0);
  vect.vector::push_back(0);
  vect.vector::push_back(0);
  vect.vector::push_back(0);

  monom mon(8,5,vect);
  mon.print_monom();

  int i, j;





//  for (i=2; i<7; ++i) {
//    for (j=1; j<i; ++j) {
//      mon.set_dg(i, j);
//      mon.print_monom();
//    }
//  }

  return 0; 
}
