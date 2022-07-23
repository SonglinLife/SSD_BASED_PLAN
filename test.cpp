#include <map>
#include <unordered_map>
#include <iostream>
#include "vector"
#include <random>
#include "queue"
struct unfill{
    unsigned pid;
    unsigned cnt;
    unfill(int pid, int cnt):pid(pid), cnt(cnt){}
    friend bool operator < (unfill a, unfill b){
        return a.cnt > b.cnt;
    }
};
int main(int argc, char const *argv[])
{
    std::unordered_map<int, int>a;
    std::priority_queue<unfill, std::vector<unfill>> p;
    a.erase(1);
    a.erase(1);
    p.push({1, 2});
    p.push({1, 3});
    std::cout<< p.top().cnt;
    std::cout << sizeof(unsigned long int) << std::endl;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    for(int i=0; i<1000000; i++){
        dis(gen);
    }
    
    return 0;
}
