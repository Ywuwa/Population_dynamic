#include "clustering.h"

#include <map>
#include <unordered_map>
#include <random>
//============================================== EUCLID DISTANCE ====================================================
double dist(const double* x, const double* y, const double xp, const double yp, const int dim) {
    double res = 0;
    for (int i = 0; i < dim; i++) {
        const double buf = x[i] - y[i];
        res += pow(buf, dim);
    }
    //const double buf = xp - yp;
    //res += pow(buf, dim+1);
    res = pow(res, 1.0/(dim)); // idk euclid or minkovsky dist
    return res;
}

//====================================== WEIGHT-CORRESPONDED NEW POPULATION INIT ====================================
void weights(const double N, const std::vector<population> &poplist, population &new_p) {
    const int dim = poplist[0].dim;
    double coef = poplist[0].N / N;

    std::unique_ptr<double []> x = std::make_unique<double []> (dim);
    for (int i = 0; i < dim; i++) {
        x[i] = poplist[0].x[i] * coef;
    }

    double p = poplist[0].p * coef;
    for (size_t i = 1; i < poplist.size(); i++) {
        coef = poplist[i].N / N;
        p += poplist[i].p * coef;
        for (int j = 0; j < dim; j++) {
            x[j] += poplist[i].x[j] * coef;
        }
    }

    new_p.init(dim, N, p, x.get());
}

//============================================== CLUSTERING =========================================================
void pop_clust(std::vector<population> &poplist, double eps) {
    size_t j = poplist.size() - 1;
    while (poplist[j].N < 1e-06) {
        poplist.pop_back();
        j -= 1;
    }

    if (poplist.size() > 1) {
        bool cflag = false;
        double N = 0.0;
        std::unordered_map<unsigned int, std::vector<population>> clust;       // dictionary with int keys and vector values

        j = poplist.size() - 1;
        clust[0].push_back(poplist[j]); // set the first pair in dict

        for (j = poplist.size() - 1; j > 0; j--) {
            // fill the dictionary (clustering)
            unsigned int key_closest = 0;
            double min_dist= 0.1;
            if (poplist[j-1].N < 1e-06) {
                poplist.pop_back();
            }
            else {
                for (auto &i: clust) {
                    const double d = dist(i.second[0].x.get(), poplist[j-1].x.get(), i.second[0].p, poplist[j-1].p, i.second[0].dim);
                    if (d < eps) {
                        //i.second.push_back(poplist[j-1]);
                        cflag = true;
                        if (d<min_dist) {min_dist = d; key_closest = i.first;}
                    }
                }
                if (cflag) {
                    clust[key_closest].push_back(poplist[j-1]);
                }
                else {
                    clust[clust.size()].push_back(poplist[j-1]);
                }
            }
            cflag = false;
            key_closest = 0;
            min_dist = 0.1;
        }

        poplist.clear(); // clear the vector in order to re-fill it
        //printf("clusters amount: %u\n", clust.size());
        for (auto &i: clust) {
            N = i.second[0].N;
            for (size_t k = 1; k < i.second.size(); k++)
                N += i.second[k].N;

            population p;
            weights(N, i.second, p);
            poplist.push_back(std::move(p));
        }
    }
}
//============================================ POPULATION RANDOM CHOICE ==========================================
int pop_choice(std::vector<population> &poplist)
{
    std::srand(time(0));
    double pop_sum = 0; // total population
    for (size_t i=0; i<poplist.size(); i++){
        pop_sum += poplist[i].N;
    }

    size_t n = poplist.size();
    std::vector<double> probab_vec(n,0);
    std::vector<double> choice_vec(n,0);
    probab_vec[0] = (poplist[0].N)/pop_sum;
    choice_vec[0] = probab_vec[0];
    for (size_t i=1; i<n; i++){
        probab_vec[i] = (poplist[i].N)/pop_sum;
        choice_vec[i] = choice_vec[i-1] + probab_vec[i];
    }

    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    const double res = distribution(generator);
    if (res<choice_vec[0]){
        //printf("res: %lf, choice_vec[0]: %lf\n", res, choice_vec[0]);
        return 0;
    }
    for (size_t i=1; i<n; i++){
        if (res > choice_vec[i-1] && res <= choice_vec[i]){
            return i;
        }
    }
    return 0;
}

//============================================ CLUSTERING 2 =========================================================
void pop_clust2(std::vector<population> &poplist, double eps) {
    size_t j = poplist.size() - 1;
    while (poplist[j].N < 1e-06) {
        poplist.pop_back();
        j -= 1;
    }

    if (poplist.size() > 1) {
        double min_dist = 0.1;
        size_t key_closest;
        bool cflag = false;
        size_t len = poplist.size();
        std::vector<population> buflist = poplist;
        std::vector<bool> flag(len, true);
        poplist.clear();
        for (j = 0; j < len; j++) {
            if (flag[j]){
            for (size_t i=j+1; i<len; i++){
                if (flag[i]){
                const double d = dist(buflist[i].x.get(), buflist[j].x.get(), buflist[i].p, buflist[j].p, buflist[j].dim);
                if (d < eps && fabs(buflist[i].p - buflist[j].p) < 0.01 && i!=j) {
                    if (d<min_dist) {min_dist = d; key_closest = i;}
                    key_closest = i;
                    cflag=true;
                    //break;
                }
                }
            }
            }
            if (cflag){
                const double del = buflist[key_closest].N + buflist[j].N;
                buflist[key_closest].p = (buflist[key_closest].p*buflist[key_closest].N + buflist[j].p*buflist[j].N)/del;
                for (int i=0; i<buflist[j].dim; i++){
                     buflist[key_closest].x[i] = (buflist[key_closest].x[i]*buflist[key_closest].N + buflist[j].x[i]*buflist[j].N)/del;
                }
                buflist[key_closest].N = del;
                flag[j]=false;
            }
            cflag = false;
            min_dist = 0.1;
        }
        for (j=0; j<len; j++){
            if (flag[j]) poplist.push_back(buflist[j]);
        }
    }
}
