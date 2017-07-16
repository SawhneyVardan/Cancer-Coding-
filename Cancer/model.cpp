//
//  main.cpp
//  Cancer
//
//  Created by Vardan Sawhney on 2016-12-04.
//  Copyright © 2016 Vardan Sawhney. All rights reserved.
//
#include <vector>
#include <algorithm>
#include <stdlib.h>

using namespace std;

struct cell {//defining cell
    int place;
    char p;
    bool is_stem;
};

static const int N = 2000; //lattice size
bool lattice[N*N] = {false}; //empty lattice
vector<cell> cells; //vector containing all cells present in the system

static const int indcNeigh[] = {-N-1, -N, -N+1, -1, 1, N-1, N, N+1};//neighborhood

char pmax=10; //proliferation capacity
double pDiv=1./24.; //division probability
double alpha=0.05; //spontaneous death probability
double ps=0.05; //probability of symmetric division
double pmig=10./24.; //probability of migration


void initialize() {
    for (int i=0; i<N; i++) {lattice[i]=true;};//filling left
    for (int i=0; i<N*N; i=i+N) {lattice[i]=true;};//filling top
    for (int i=N-1; i<N*N; i=i+N) {lattice[i]=true;};//filling bottom
    for (int i=N*(N-1); i<N*N; i++) {lattice[i]=true;};//filling right
    
    lattice[N/2*N+N/2] = true; //initial cell in the middle
    cell initialCell = {N/2*N+N/2,pmax,true};//initial cell definition
    cells.push_back(initialCell);
}

int returnEmptyPlace(int indx) {
    int neigh[8], nF = 0;
    for(int j=0;j<8;j++) {//searching through neighborhood
        if (!lattice[indx+indcNeigh[j]]) {//if free spot
            neigh[nF] = indx+indcNeigh[j]; //save the index
            nF++; //increase the number of found free spots
        }
    }
    if(nF) {//selecting free spot at random
        return neigh[rand() % nF];
    } else {//no free spot
        return 0;
    }
}

void simulate(int nSteps) {
    
    vector<cell> cellsTmp;
    int newSite;
    cell currCell, newCell;
    
    for (int i=0; i<nSteps; i++) {
        random_shuffle(cells.begin(), cells.end()); //shuffling cells
        while (!cells.empty()) {
            currCell=cells.back(); //pick the cell
            cells.pop_back();
            newSite = returnEmptyPlace(currCell.place);
            
            if (newSite) {//if there is a new spot
                newCell = currCell;
                newCell.place = newSite;
                if ((double)rand()/(double)RAND_MAX < pDiv) {
                    if (currCell.is_stem) {
                        lattice[newSite]=true;
                        cellsTmp.push_back(currCell);
                        if ((double)rand()/(double)RAND_MAX > ps) {//asymmetric division
                            newCell.is_stem = false;
                        }
                        cellsTmp.push_back(newCell);
                    } else if (currCell.p > 0 && (double)rand()/(double)RAND_MAX > alpha) {
                        currCell.p--;
                        newCell.p--;
                        lattice[newSite] = true;
                        cellsTmp.push_back(currCell);
                        cellsTmp.push_back(newCell);
                    } else {
                        lattice[currCell.place] = false;
                    }
                } else if ((double)rand()/(double)RAND_MAX < pmig) {
                    lattice[currCell.place] = false;
                    lattice[newSite] = true;
                    cellsTmp.push_back(newCell);
                } else {//doing nothing
                    cellsTmp.push_back(currCell);
                }
            } else {//no free spot
                cellsTmp.push_back(currCell);
            }
        }
        cells.swap(cellsTmp);
    }
}

int main() {
    srand(time_t(NULL)); //initialize random number generator
    initialize(); //initialize CA
    simulate(24*30*6);
    return 0;
}
