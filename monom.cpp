#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// zzbar /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

class zzbar {
  public:
    zzbar();
    zzbar(int d, int n);
    zzbar(const zzbar& copy);
    void set_dg(int d, int n);
    int  alpha();
    int  beta();
    int  prefix();
    int  get_deg();
    int  get_num();
    void print_zzbar();
    bool operator<(zzbar ob2);
    bool operator==(zzbar ob2);
  protected:
    int deg;  // degree of monom z^{\alpha}\bar{z}^{\beta}
    int num;  // numer of monom among all monoms degree deg
};


zzbar::zzbar() {
  deg = 2;
  num = 1;
}

zzbar::zzbar(int d, int n) {
  if (d > 1 && n < d) {
    deg = d;
    num = n;
  }
}

zzbar::zzbar(const zzbar& copy) {
  deg = copy.deg;
  num = copy.num;
}

// set deg and num
void zzbar::set_dg(int d, int n) {
  deg = d;
  num = n;
}

// Function finding \alpha
int zzbar::alpha() {
  int halfnum = (1 + num) / 2;
  int alp = deg - halfnum;
  return alp;
}

// Function finding \beta
int zzbar::beta() {
  int bet = (1 + num) / 2;
  return bet;
}

// Function finding prefix
int zzbar::prefix() {
  int prefix;
  if (deg%2 == 0 && num == deg-1) 
          prefix = 3;    // without prefix
  else {
    if(num%2 == 1) 
          prefix = 1;    // prefix = Re
    if(num%2 == 0) 
          prefix = 2;    // prefix = Im
  }
  return  prefix;
}

int zzbar::get_deg() {
  return deg;
}

int zzbar::get_num() {
  return num;
}

// Operator < comparing two monoms zzbar
bool zzbar::operator<(zzbar ob2) {
  bool answ = true;
  if (deg < ob2.deg)
    answ = true;
  else if (deg == ob2.deg && num > ob2.num)
    answ = true;
  else 
    answ = false;

  return answ;
}

// Operator == comparing two monoms zzbar
bool zzbar::operator==(zzbar ob2) {
  bool answ = true;
  if (deg == ob2.deg && num == ob2.num) {
    answ = true;
  }
  else answ = false;

  return answ;
}

// Function printing monom zzbar
void zzbar::print_zzbar() {
  int pr  = prefix();
  int alp = alpha();
  int bet = beta();

  if (pr == 1) 
    std::cout << "2\\mathrm{Re}\\left(";
  else if (pr == 2) 
    std::cout << "2\\mathrm{Im}\\left(";

  if (alp == 1) 
    std::cout << "z";
  else 
    std::cout << "z^{" << alp << "}";

  if (bet == 1) 
    std::cout << "\\bar{z}";
  else 
    std::cout << "\\bar{z}^{" << bet << "}";

  if(pr < 3) 
    std::cout << "\\right)";
}


//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// monom /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

class monom : public zzbar {
  public:
    monom();
    monom(int d, int n, const std::vector<int> &vec);
    monom(int d, int n, const std::vector<int> &vec, const std::vector<int> &deg_wei);
    monom(const monom& copy);
    void set(int d, int n, const std::vector<int> &vec);
    void print_monom();
    void set_index(int ind);
    void set_weight(int wei);                
    void regul();
    bool operator<(const monom &ob2) const;
    bool operator==(const monom &ob2) const;
    int  get_index();
    int  get_weight();
    int  get_deguwei(const int i);
    void get_deguwei(std::vector<int> &vec);
    void get_degu(std::vector<int> &vec);
    int  get_degusize();

  private:
    int index;
    int weight;
    std::vector<int> degu;      // degu = <g_1, ..., g_s>, where monom = zzbar * u_1^{g_1} * ... * u_s^{g_s}
    std::vector<int> degu_wei;  // degu_wei = <w_0, w_1, ..., w_s>, where w_i = weight of cutting monom = zzbar * u_1^{g_1} * ... * u_s^{g_i}

    static bool less_vector(const std::vector<int> &v1, const std::vector<int> &v2);
};


monom::monom() {
  deg    = 2;
  num    = 1;
  index  = 1;
  weight = 2;
  degu_wei.push_back(2);
}

monom::monom(int d, int n, const std::vector<int> &vec)
    : zzbar(d, n) {
  degu = vec;
  while (!degu.empty() && degu.back() == 0) {
    degu.pop_back();
  }
  index  = 1;
  weight = 2;
}

monom::monom(int d, int n, const std::vector<int> &vec, const std::vector<int> &deg_wei)
    : zzbar(d, n) {
  degu = vec;
  while (!degu.empty() && degu.back() == 0) {
    degu.pop_back();
  }
  degu_wei = deg_wei;
  index    = 1;
  weight   = 2;
}

monom::monom(const monom& copy) {
  deg      = copy.deg;
  num      = copy.num;
  index    = copy.index;
  weight   = copy.weight;
  degu     = copy.degu;
  degu_wei = copy.degu_wei;
}

void monom::set(int d, int n, const std::vector<int> &vec) {
  degu = vec;
  while (!degu.empty() && degu.back() == 0) {
    degu.pop_back();
  }
  index  = 1;
  weight = 2;
  deg = d;
  num = n;
}

void monom::set_index(int ind) {
  index = ind;
}

void monom::set_weight(int wei) {
  weight = wei;
}

bool monom::less_vector(const std::vector<int> &v1, const std::vector<int> &v2) {
  bool answ = true;
  int i = 0;
  int flag = 0;
  int n = std::min(v1.size(), v2.size());

  while (flag == 0 && i<n) {
    if (v1[i] < v2[i]) flag = 1;
    ++i;
  }

  if (flag == 1) answ = true;
  else           answ = false;

  return answ;
}

void monom::regul() {
  while (!degu.empty() && degu.back() == 0) {
    degu.pop_back();
  }
  int n = degu.size();
  int m = degu_wei.size();
  if(m > n+1) {
    int k = m-n-1;
    for(int i = 0; i < k; ++i) {
      degu_wei.pop_back();
    }
  }
}

// Operator < comparing two monoms
bool monom::operator<(const monom &ob2) const {
  bool answ = true;

  if (deg < ob2.deg) 
    answ = true;
  else if (deg == ob2.deg && num > ob2.num) 
    answ = true;
  else if (deg == ob2.deg && num == ob2.num && less_vector(degu, ob2.degu)) 
    answ = true;
  else 
    answ = false;

  return answ;
}

// Operator == comparing two monoms
bool monom::operator==(const monom &ob2) const {
  bool answ = true;
  std::vector<int> vec_l = degu;
  std::vector<int> vec_r = ob2.degu;

  while (!vec_l.empty() && vec_l.back() == 0) {
    vec_l.pop_back();
  }
  while (!vec_r.empty() && vec_r.back() == 0) {
    vec_r.pop_back();
  }

  if (deg == ob2.deg && num == ob2.num && vec_r == vec_l) 
    answ = true;
  else 
    answ = false;

  return answ;
}

// Function printing monom
void monom::print_monom() {
  int n = degu.size();

  std::cout << "$$v_{" << index << "} = ";
  print_zzbar();
  for (int i = 0; i < n; ++i) {
    if (degu[i] > 1) 
      std::cout << "u_{" << i+1 << "}^{" << degu[i] << "}";
    else if (degu[i] == 1) 
      std::cout << "u_{" << i+1 << "}";
  }
  std::cout << "$$";
  std::cout << "\n";
}

int monom::get_index() {
  return index;
}

int monom::get_weight() {
  return weight;
}

int monom::get_deguwei(int i) {
  if( i < degu_wei.size() && i >= 0) {
    return degu_wei[i];
  }
  else {
    std::cout << "error. monom::get_deguwei. i >= degu_wei.size()";
    return 0;
  }
}

void monom::get_deguwei(std::vector<int> &vec) {
  vec = degu_wei;
}

void monom::get_degu(std::vector<int> &vec) {
  vec = degu;
}

int monom::get_degusize() {
  return degu.size();
}


//////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// H_n ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

class H_n {
  public:
    H_n();
    H_n(int n, const std::vector<int> &dim);
    int get_size();
    void print_space();
    void print_space(int i);
    monom get_space(int i);
    std::vector<monom> get_space();

  private:
    int period_num;                  // number of period \in {2, 3, ...}
    std::vector<monom> space;        // all monoms of weight period_num
    std::vector<int>   dimensions;   // dimensions of spaces CH_2, ... CH_{n-1}

    static int smartloop_inc (const std::vector<int> &weight, std::vector<int> &vec, int summa);
    void spacecreator(int n, const std::vector<int> &dim);
};


H_n::H_n() {
  std::vector<int> dim;
  spacecreator(2, dim);
}

H_n::H_n(int n, const std::vector<int> &dim) {
  spacecreator(n, dim);
}

int H_n::get_size() {
  return space.size();
}

void H_n::print_space() {
  for (int i = 0; i < space.size(); ++i) {
    space[i].print_monom();
  }
}

monom H_n::get_space(int i) {
  if (space.size() > i)
    return space[i];
}

std::vector<monom> H_n::get_space() {
  return space;
}

void H_n::print_space(int i) {
  if (i < space.size())
    space[i].print_monom();
}

void H_n::spacecreator(int n, const std::vector<int> &dim) {
  if( n-2 != dim.size() ) 
    std::cout << "error. n-2 != dim.size()";
  
  int i    = 0;
  int flag = 0;
  int iwei = 0;
  int lim  = (n-4) / 2;

  period_num = n;
  dimensions = dim;

  int sumdim = 0;
  for(i=0; i<=lim; ++i) {
    sumdim += dimensions[i];
  }

  std::vector<int> weights;
  for(i = 0; i < sumdim; ++i) {
    flag = 0;
    iwei = 0;
    while (i+1 > flag) {
      flag += dimensions[iwei];
      ++iwei;
    }
    weights.push_back(iwei+1);
  }
//======== creation of vector space ========
  std::vector<int> loopv(sumdim, 0);
  for(int ideg = period_num; ideg > 1; --ideg) {
    int diff = period_num - ideg;
    for(int inum = 1; inum < ideg; ++inum) {
      loopv.assign(sumdim, 0);
      int fl = 0;
      if (diff > 0) {
        fl = smartloop_inc(weights, loopv, diff);
      } 
      while (fl == 0) {
        int sum = ideg;
        std::vector<int> deg_wei;
        deg_wei.push_back(sum);
        for (int ive = 0; ive < sumdim; ++ive) {
          sum += (loopv[ive] * weights[ive]);
          deg_wei.push_back(sum);
        }
        monom temp(ideg, inum, loopv, deg_wei);
        temp.set_weight(n);      
        temp.regul();          
        space.push_back(temp);
        fl = smartloop_inc(weights, loopv, diff);
      }
    }
  }
//==========================================
}

int H_n::smartloop_inc (const std::vector<int> &weight, std::vector<int> &vec, int summa) {
  if (weight.size() != vec.size()) 
    std::cout << "error. weight.size != vec.size";

  int end = 0;
  int n   = weight.size();
 
  if(n > 0) {
    int fl  = 0;
    vec[0] += 1;

    while ( fl == 0 ) {
      int sum = 0;
      for(int i = 0; i < n; ++i) {
        sum += weight[i] * vec[i];
      }

      if(sum < summa) 
        vec[0] += 1;
      else if (sum == summa)
        fl = 1;
      else {
        int fl2  = 0;
        int ivec = 0;
        while (fl2 == 0 && ivec < n) {
          vec[ivec]==0 ? ++ivec : fl2=1;
        }
        if( ivec == n || ivec == n-1 ) {
          fl  = 1;
          end = 1;
        }
        else {
          vec[ivec] = 0;
          vec[1+ivec] += 1;
        }
      }

    }
  }
  
  if(n == 0) 
    end = 1;

  return end;
}

//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// CH_n /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

class CH_n {
  public:
    CH_n();
    CH_n(int n);
    void print_megaspace();
    void print_info();

  private:
    std::vector< std::vector<monom> > megaspace;
    std::vector<int> dim;
    
    int  monom_index(const monom &mon);
    bool monom_check(const monom &mon);
};

CH_n::CH_n() {
  H_n H;
  megaspace.push_back( H.get_space() );
}

CH_n::CH_n(int n) {
  if(n < 2 ) 
    std::cout << "error. n < 2";

  int ind = 0;

  for(int i = 2; i < n+1; ++i) {
    H_n H(i, dim);
    dim.push_back(0);
    std::vector<monom> sp;
    for (int j = 0; j < H.get_size(); ++j) {
      if(  monom_check( H.get_space(j) )  )  {
        ind++;
        dim[i-2] += 1;
        sp.push_back( H.get_space(j) );
        sp.back().set_index(ind);
      }
    }
    megaspace.push_back( sp );
  }
}

void CH_n::print_megaspace() {
  for(int i = 0; i < megaspace.size(); ++i) {
    std::cout << "$$========== period \\ " << i+2 << " ==========$$ \n";
    for(int j = 0; j < megaspace[i].size(); ++j) {
      megaspace[i][j].print_monom();
    }
  }
}

void CH_n::print_info() {
  std::cout << "\n";
  for(int i = 0; i < megaspace.size(); ++i) {
    std::cout << "$$\\dim(CH_{" << i+2 << "}) = " << megaspace[i].size() << "$$ \n";
  }
  std::cout << "\n";
  for(int i = 0; i < megaspace.size(); ++i) {
    int max  = 0;
    int jmax = 0;
    for(int j = 0; j < megaspace[i].size(); j++) {
      if( megaspace[i][j].get_degusize() > max ) {
        max  = megaspace[i][j].get_degusize();
        jmax = j;
      }
    }
    std::cout << "$$\\dim(U, \\ CH_{" << i+2 << "}) = " << max << ", \\ \\ \\text{for example}:$$ \n";
    megaspace[i][jmax].print_monom();
  }
  std::cout << "\n";
}

int CH_n::monom_index(const monom &mon) {
  int ind  = 0;

  monom mono(mon);
  std::vector<int> deg_u;
  mono.get_degu(deg_u);

  int wei  = mono.get_deguwei( deg_u.size() ) - 2;
  
  int fl = 0;
  int i  = 0;
 
  while( fl == 0 && i < megaspace[wei].size() ) {
    if(megaspace[wei][i] == mon) {
      ind = megaspace[wei][i].get_index();
      fl = 1;
    } 
    ++i;
  }

  return  ind;
}

bool CH_n::monom_check(const monom &mon) {
  int answ = true;
  int fl   = 0;
  monom mono(mon);
  std::vector<int> degu_w;
  mono.get_deguwei(degu_w);
  std::vector<int> deg_u;
  mono.get_degu(deg_u);

  if( 2*mono.get_deg() <= mono.get_weight() ) {
    int i = 0;
    int detailed_check = 0;      
    std::vector<int> vec;
    while (fl == 0 && i < deg_u.size() ) {
      monom temp_mon(mono.get_deg(), mono.get_num(), vec, degu_w);
      vec.push_back(deg_u[i]);
      if(deg_u[i] > 0 && monom_index(temp_mon) <= i+1) {
        fl = 1;
      }
      ++i;
    } 
  }

  if(fl == 0)
    answ = true;
  else
    answ = false;

  return answ;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


int main()
{

  time_t sec1 = time(NULL);

  CH_n ch6(8);

  time_t sec2 = time(NULL);
  time_t dsec = sec2 - sec1;

  ch6.print_info();
  ch6.print_megaspace();

  std::cout << "Elapsed time: " << dsec << " sec" << std::endl;

  return 0;
}
