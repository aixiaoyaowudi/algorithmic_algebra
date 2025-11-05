#include <polynomials.h>
#include <classes.h>
#include <random>
#include <chrono>
std::mt19937 rnd(std::chrono::high_resolution_clock::now().time_since_epoch().count());

int main()
{
    constexpr int N(30);
    std::vector<Z> v(N-1);
    int g(std::uniform_int_distribution<int>(N,N*N)(rnd));
    for(auto &t:v) if(std::uniform_int_distribution<int>(0,1)(rnd))
    {
        t.set(g);
        for(int _(0);_<5;++_)
        {
            t = t * std::uniform_int_distribution<int>(1,N)(rnd);
        }
        if(std::uniform_int_distribution<int>(0,1)(rnd)) t = - t;
    }
    // v[0] = 6 *g; v[1] = 15*g;v[2] = 10*g;
    // v[0] = 3 * g; v[1] = 5 *g;
    // v.emplace_back();
    std::cout<<g<<std::endl;
    auto [res, _] = detach(v, g * g);
    std::cout<<_<<std::endl;
    for(int i(0); i<N;++i)
    {
        std::cout<<"("<<res[i]<<")*("<<v[i]<<")"<<("+="[i==N-1]);
    }
    std::cout<<std::endl;
    // std::cout<<g * g<<std::endl;
    return 0;
}