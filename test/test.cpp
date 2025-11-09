#include <polynomials.h>
#include <classes.h>
#include <random>
#include <chrono>
std::mt19937 rnd(std::chrono::high_resolution_clock::now().time_since_epoch().count());

int main()
{
    // constexpr int N(30);
    // std::vector<Z> v(N-1);
    // int g(std::uniform_int_distribution<int>(N,N*N)(rnd));
    // for(auto &t:v) if(std::uniform_int_distribution<int>(0,1)(rnd))
    // {
    //     t.set(g);
    //     for(int _(0);_<5;++_)
    //     {
    //         t = t * std::uniform_int_distribution<int>(1,N)(rnd);
    //     }
    //     if(std::uniform_int_distribution<int>(0,1)(rnd)) t = - t;
    // }
    // // v[0] = 6 *g; v[1] = 15*g;v[2] = 10*g;
    // // v[0] = 3 * g; v[1] = 5 *g;
    // // v.emplace_back();
    // std::cout<<g<<std::endl;
    // auto [res, _] = detach(v, g * g);
    // std::cout<<_<<std::endl;
    // for(int i(0); i<N;++i)
    // {
    //     std::cout<<"("<<res[i]<<")*("<<v[i]<<")"<<("+="[i==N-1]);
    // }
    // std::cout<<std::endl;
    // std::cout<<g * g<<std::endl;
    // using poly = multivariate_polynomial_ring<Q, term_ordering::lexicographic_ordering>;
    // ideal_solver<Z> ih({Z(1)});
    // auto [A,B]=ih.detach(Z(1));
    // std::cerr<<B<<std::endl;
    using poly = multivariate_polynomial_ring<Z, term_ordering::lexicographic_ordering>;
    // poly a,b;
    poly a,b,c;
    // a.insert({std::make_pair(Hterm({0,0,1,0,0}),Q(1)),
    //           std::make_pair(Hterm({1,0,0,0,0}),Q(-1)),
    //           std::make_pair(Hterm({0,1,0,0,0}),Q(-1))});
    // b.insert({std::make_pair(Hterm({0,0,0,1,0}),Q(1)),
    //           std::make_pair(Hterm({2,0,0,0,0}),Q(-1)),
    //           std::make_pair(Hterm({1,1,0,0,0}),Q(-2))});
    // c.insert({std::make_pair(Hterm({0,0,0,0,1}),Q(1)),
    //           std::make_pair(Hterm({3,0,0,0,0}),Q(-1)),
    //           std::make_pair(Hterm({2,1,0,0,0}),Q(-3))});
    a.insert({std::make_pair(Hterm({0,0,1,0,0}),Z(1)),
              std::make_pair(Hterm({1,0,0,0,0}),Z(-1)),
              std::make_pair(Hterm({0,1,0,0,0}),Z(-1))});
    b.insert({std::make_pair(Hterm({0,0,0,1,0}),Z(1)),
              std::make_pair(Hterm({2,0,0,0,0}),Z(-1)),
              std::make_pair(Hterm({1,1,0,0,0}),Z(-2))});
    c.insert({std::make_pair(Hterm({0,0,0,0,1}),Z(1)),
              std::make_pair(Hterm({3,0,0,0,0}),Z(-1)),
              std::make_pair(Hterm({2,1,0,0,0}),Z(-3))});
    // std::cout << term_ordering::lexicographic_ordering()(Hterm({1}),Hterm({0})) << " " << term_ordering::lexicographic_ordering()(Hterm({0}),Hterm({1}))<<std::endl;
    // // std::cout << a << " " << b <<std::endl;
    // // a.insert(std::make_pair(Hterm({1}), Q(1)));
    // // std::cout << a << " " << b <<std::endl;
    // // a.insert(std::make_pair(Hterm({0}), Q(1)));
    // // std::cout << a << " " << b <<std::endl;
    // // b.insert(std::make_pair(Hterm({2}), Q(1)));
    // // std::cout << a << " " << b <<std::endl;
    // // b.insert(std::make_pair(Hterm({0}), Q(1)));
    // // std::cout << a << " " << b <<std::endl;
    // a.insert({std::make_pair(Hterm({2}),Z(1)),std::make_pair(Hterm({0,1}),Z(-1))});
    // b.insert({std::make_pair(Hterm({1,1}),Z(1)),std::make_pair(Hterm({0}),Z(-1))});
    std::cout << a << " " << b << " " << c << std::endl;
    // std::cout << a << " " << b << std::endl;
    auto T = grobner_basis(std::vector<poly>({a,b,c}));
    // auto T = grobner_basis(std::vector<poly>({a,b}));
    int idx(0);
    for(auto [g, coef] : T)
    {
        std::cout << "Basis #" << (++idx) << " " << g << " ";
        for(auto v : coef) std::cout << v << " ";
        std::cout<<std::endl;
    }
    // std::vector<Z> c = {8, 18, 12, 0};
    // syzygy_solver<Z> ch(c);
    // for(auto val:ch.basis)
    // {
    //     for(auto v:val) std::cerr<<v<<",";
    //     std::cerr<<std::endl;
    // }
    return 0;
}