#include "Haplotype.hpp"

#define SAHAP_CHROMOSOME_DEBUG_MEC 0 // Not sure what this is meant to debug.... chromosome?
#if SAHAP_CHROMOSOME_DEBUG_MEC
#include <csignal>
#endif
#define SAHAP_CHROMOSOME_ALT_MEC 0
#define SMALL_ENOUGH_TO_IGNORE 1e-10

namespace SAHap {

    Haplotype::Haplotype(dnapos_t numSites, unsigned ploidyCount)
        : numSites(numSites), total_mec(0), window_mec(0), isitecost(0)
    {
        this->ploidyCount = ploidyCount;
        this->solution = vector<int>(this->numSites, -1);
        this->weights = vector<vector<int>>(this->numSites);
        this->siteCoverages = vector<dnacnt_t>(this->numSites);
        this->window.start = 0;
        this->window.end = this->numSites;

        for (dnapos_t i = 0; i < this->numSites; ++i) {
            this->solution[i] = -1;
            this->weights[i] = vector<int>(this->ploidyCount, 0);
            this->siteCoverages[i] = 0;
        }
    }
