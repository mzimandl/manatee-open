//  Copyright (c) 2016-2023  Milos Jakubicek

#include <cmath>
#include <finlib/regpref.hh>
#include "subcorp.hh"
#include "keyword.hh"

bool check_string (string &str, vector<regexp_pattern*> pos_filters, vector<regexp_pattern*> neg_filters,
                   unordered_set<string> &blist, unordered_set<string> &wlist)
{
    for (auto it = neg_filters.begin(); it != neg_filters.end(); it++) {
        if ((*it)->match(str.c_str()))
            return false;
    }
    for (auto it = pos_filters.begin(); it != pos_filters.end(); it++) {
        if (!(*it)->match(str.c_str()))
            return false;
    }
    if (!blist.empty() && blist.find(str) != blist.end())
        return false;
    if (!wlist.empty() && wlist.find(str) == wlist.end())
        return false;
    return true;
}

class AllowMissingFrequency : public Frequency {
    unique_ptr<Frequency> src;
public:
    AllowMissingFrequency (WordList *wl, const char *frqtype) {
        try {
            src.reset(wl->get_stat(frqtype));
        } catch (FileAccessError&) {}
    }
    double freq(int id) {return src ? src->freq(id) : -1;}
};

Keyword::Keyword (Corpus *c1, Corpus *c2, WordList *wl1, WordList *wl2, float N,
             unsigned maxlen, int minfreq, int maxfreq, unordered_set<string> &blacklist,
             unordered_set<string> &whitelist, const char *frqtype, vector<string> &addfreqs,
             vector<string> pos_regex_filters, vector<string> neg_regex_filters, FILE* progress)
            : curr (0), totalcount(0), totalfreq1(0), totalfreq2(0)
{
    auto get_data_index = [](const string& scoretype, const int addfreqs_size) -> int {
        if (scoretype == "logL")
            return 2 * addfreqs_size + 5;
        else if (scoretype == "chi2")
            return 2 * addfreqs_size + 6;
        else if (scoretype == "din")
            return 2 * addfreqs_size + 7;
        else
            return 2 * addfreqs_size + 4;
    };

    string str_params(frqtype), ftype, sortby, filterby, substr_params;
    double filter_min = std::nan(""), filter_max = std::nan("");
    int n = str_params.find(";");
    if (n == -1) {
        ftype = str_params;
    } else {
        ftype = str_params.substr(0, n);
        substr_params = str_params.substr(n+1);

        int m = substr_params.find(";");
        if (m == -1) {
            sortby = substr_params;
        } else {
            sortby = substr_params.substr(0, m);
            substr_params = substr_params.substr(m+1);
            if (!substr_params.empty()) {
                size_t slash1 = substr_params.find("/");
                size_t slash2 = substr_params.rfind("/");
                if (slash1 != string::npos && slash2 != string::npos && slash1 != slash2) {
                    filterby = substr_params.substr(0, slash1);
                    string substr_min = substr_params.substr(slash1 + 1, slash2 - slash1 - 1);
                    if (!substr_min.empty())
                        filter_min = stod(substr_min);
                    string substr_max = substr_params.substr(slash2 + 1);
                    if (!substr_max.empty())
                        filter_max = stod(substr_max);
                } else {
                    throw new CorpInfoNotFound("Invalid filterby format for Keyword");
                }
            }
        }
    }

    const int sort_index = get_data_index(sortby, addfreqs.size());
    const int filter_index = get_data_index(filterby, addfreqs.size());

    vector<regexp_pattern*> pos_regpats;
    for (auto it = pos_regex_filters.begin(); it != pos_regex_filters.end(); it++) {
        if (!(*it).size())
            continue;
        regexp_pattern *re = new regexp_pattern((*it).c_str(), wl1->locale, wl1->encoding);
        if (re->compile())
            throw new CorpInfoNotFound("Invalid regular expression for Keyword");
        pos_regpats.push_back(re);
    }

    vector<regexp_pattern*> neg_regpats;
    for (auto it = neg_regex_filters.begin(); it != neg_regex_filters.end(); it++) {
        if (!(*it).size())
            continue;
        regexp_pattern *re = new regexp_pattern((*it).c_str(), wl1->locale, wl1->encoding);
        if (re->compile())
            throw new CorpInfoNotFound("Invalid regular expression for Keyword");
        neg_regpats.push_back(re);
    }

    strid_gen *it = NULL;
    if (c1->get_confpath() == c2->get_confpath())
        it = new strid_gen (wl1->dump_str());
    else
        it = new strid_gen (new WordListLeftJoin (wl1->attr_path, wl2->attr_path));

    const float c1size = c1->search_size();
    const float c2size = c2->search_size();
    const float c12size = c1size + c2size;
    const float wl1_maxid = wl1->id_range();
    Frequency *stat1 = wl1->get_stat(ftype.c_str());
    Frequency *stat2 = wl2->get_stat(ftype.c_str());
    vector<AllowMissingFrequency> addfreqs1, addfreqs2;
    for (auto it = addfreqs.begin(); it != addfreqs.end(); it++) {
        addfreqs1.emplace_back(wl1, (*it).c_str());
        addfreqs2.emplace_back(wl2, (*it).c_str());
    }

    while (!it->end()) {
        string str;
        int id1, id2;
        it->next (str, id1, id2);
        if (progress && id1 % 100000 == 0)
            fprintf(progress, "\r%.2f %%", id1 / wl1_maxid * 100);
        NumOfPos f1 = stat1->freq(id1);
        NumOfPos f2;
        if (id2 == -1)
            f2 = 0;
        else
            f2 = stat2->freq(id2);
        if (minfreq && f1 < minfreq)
            continue;
        if (maxfreq && f1 > maxfreq)
            continue;
        totalcount++;
        totalfreq1 += f1;
        totalfreq2 += f2;
        float fpm1 = f1 * 1000000 / c1size;
        float fpm2 = f2 * 1000000 / c2size;
        // adding 4 additional fields to `freqs`
        double *freqs = new double[2*addfreqs.size() + 4 + 4];
        freqs[0] = f1; freqs[1] = f2;
        freqs[2] = fpm1; freqs[3] = fpm2;
        for (unsigned i = 0; i < addfreqs.size(); i++) {
            freqs[2*i+4] = addfreqs1[i].freq(id1);
            freqs[2*i+5] = id2 == -1 ? 0 : addfreqs2[i].freq(id2);
        }

        double e1 = c1size * (f1 + f2) / c12size;
        double e2 = c2size * (f1 + f2) / c12size;
        freqs[2*addfreqs.size() + 4] = (fpm1 + N) / (fpm2 + N); // manatee score
        freqs[2*addfreqs.size() + 5] = 2 * (f1*log(f1/e1) + f2*log(f2/e2)); // logL score
        freqs[2*addfreqs.size() + 6] = (f1-e1)*(f1-e1)/e1 + (f2-e2)*(f2-e2)/e2; // chi2 score
        freqs[2*addfreqs.size() + 7] = 100 * ((fpm1 - fpm2) / (fpm1 + fpm2)); // DIN size effect

        float score;
        if (isnan(freqs[sort_index])) {
            score = 0;
        } else {
            score = freqs[sort_index];
        }

        if (!isnan(filter_min) && (isnan(freqs[filter_index]) || freqs[filter_index] < filter_min)) {
            delete[] freqs;
            totalcount--;
            totalfreq1 -= f1;
            totalfreq2 -= f2;
            continue;
        }
        if (!isnan(filter_max) && (isnan(freqs[filter_index]) || freqs[filter_index] > filter_max)) {
            delete[] freqs;
            totalcount--;
            totalfreq1 -= f1;
            totalfreq2 -= f2;
            continue;
        }

        if (heap.size() < maxlen) {
            if (!check_string (str, pos_regpats, neg_regpats, blacklist, whitelist)) {
                delete[] freqs;
                totalcount--;
                totalfreq1 -= f1;
                totalfreq2 -= f2;
                continue;
            }
            kwitem *k = new kwitem (id1, id2, score, str, freqs);
            heap.push_back (k);
            push_heap (heap.begin(), heap.end(), kwitem_cmp());
        } else if (heap.front()->score < score) {
            if (!check_string (str, pos_regpats, neg_regpats, blacklist, whitelist)) {
                delete[] freqs;
                totalcount--;
                totalfreq1 -= f1;
                totalfreq2 -= f2;
                continue;
            }
            pop_heap (heap.begin(), heap.end(), kwitem_cmp());
            delete heap.back();
            heap.pop_back();
            kwitem *k = new kwitem (id1, id2, score, str, freqs);
            heap.push_back (k);
            push_heap (heap.begin(), heap.end(), kwitem_cmp());
        } else
            delete[] freqs;
    }
    delete it;
    delete stat1; delete stat2;
    for (auto it = pos_regpats.begin(); it != pos_regpats.end(); it++)
        delete *it;
    for (auto it = neg_regpats.begin(); it != neg_regpats.end(); it++)
        delete *it;
    if (progress)
        fprintf(progress, "\r100 %%     \n");
    sort_heap (heap.begin(), heap.end(), kwitem_cmp());
}

// vim: ts=4 sw=4 sta et sts=4 si cindent tw=80:
