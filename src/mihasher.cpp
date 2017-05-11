#include <algorithm>
#include "mihasher.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <inttuple.h>
#include <string>
using namespace std;

/*
 * Inputs: query, numq, dim1queries
 */

/*
 * Outputs: results, numres, stats
 *
 *   results: an array of indices (1-based) of the K-nearest neighbors
 *   for each query. So the array includes K*numq uint32 integers.
 *
 *   numres: includes the number of database entries that fall at any
 *   specific Hamming distance from the query until the K nearest
 *   neighbors are reached. So from this array you can figure out the
 *   Hamming distances of the K-nearest neighbors.
 */

//extern int num_comparisons;
bool comp(radius R,radius Rp)
{
  //num_comparisons++;
  double sim1 = (double)(R.querynorm-R.r1)/(sqrt(R.querynorm-R.r1+R.r2));
  double sim2 = (double)(Rp.querynorm-Rp.r1)/(sqrt(Rp.querynorm-Rp.r1+Rp.r2));
  // cout<<"Querynorm = "<< R.querynorm <<" R.r1 = "<<R.r1 <<" R.r2 =" << R.r2 <<" sim1 = :"<<sim1<<" sim2 = :"<<sim2<<"\n";
  // cout<<"Querynorm = "<< Rp.querynorm <<" Rp.r1 = "<<Rp.r1 <<" Rp.r2 =" << Rp.r2 <<" sim1 = :"<<sim1<<" sim2 = :"<<sim2<<"\n";
	return sim1<sim2;
}

void mihasher::batchquery(UINT32 *results, UINT32 *numres, qstat *stats, UINT8 *queries, UINT32 numq, int dim1queries, int& max)
{
	counter = new bitarray;
	counter->init(N);

	UINT32 *res  = new UINT32[K*(D+1)];
	UINT64 *chunks = new UINT64[m];
	max = 0;
	UINT32 *presults = results;
	UINT32 *pnumres = numres;
	UINT64 num_comp = 0,num_buck_check = 0;
	qstat *pstats = stats;
	int maxheap = 0;
	bool secondcase = false;
	int num_sec_case = 0;
	double rati;
	UINT8* 	pq = queries;

	for (int i=0; i<numq; i++) {
	  query(presults, pnumres, pstats, pq, chunks, res);
		presults += K;
		max = maxheap + max;
		pnumres += B+1;
		pstats ++;
		pq += dim1queries;
		
	}
	delete [] res;
	delete [] chunks;

	delete counter;
}


// Temp variables: chunks, res -- I did not want to malloc inside
// query, so these arrays are passed from outside

void mihasher::query(UINT32 *results, UINT32* numres, qstat *stats, UINT8 *query, UINT64 *chunks, UINT32 *res)
{	
	UINT32 maxres = K ? K : N;			// if K == 0 that means we want everything to be processed.
	// So maxres = N in that case. Otherwise K limits the results processed.

	UINT32 n = 0; 				// number of results so far obtained (up to a distance of s per chunk)
	UINT32 nc = 0;				// number of candidates tested with full codes (not counting duplicates)
	UINT32 nd = 0;                      	// counting everything retrieved (duplicates are counted multiple times)
	UINT32 nl = 0;				// number of lookups (and xors)
	UINT32 *arr;
	int size = 0;
	UINT32 index;
	int hammd;
	clock_t start, end;

	start = clock();

	counter->erase();
	memset(numres, 0, (B+1)*sizeof(*numres));

	split(chunks, query, m, mplus, b);

	int s;			// the growing search radius per substring

	int curb = b;		// current b: for the first mplus substrings it is b, for the rest it is (b-1)


	for (s = 0; s <= d && n < maxres; s++) {
		for (int k=0; k<m; k++) {
			if (k < mplus)
				curb = b;
			else
				curb = b-1;
			UINT64 chunksk = chunks[k];
			nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s

			UINT64 bitstr = 0; 			// the bit-string with s number of 1s
			for (int i=0; i<s; i++)
				power[i] = i;			// power[i] stores the location of the i'th 1
			power[s] = curb+1;			// used for stopping criterion (location of (s+1)th 1)

			int bit = s-1;			// bit determines the 1 that should be moving to the left
			// we start from the left-most 1, and move it to the left until it touches another one
			//cout<<"here2\n";
			while (true) {			// the loop for changing bitstr
				if (bit != -1) {
					bitstr ^= (power[bit] == bit) ? (UINT64)1 << power[bit] : (UINT64)3 << (power[bit]-1);
					power[bit]++;
					bit--;
				} else {
					/* the binary code bitstr is available for processing */
					arr = H[k].query(chunksk ^ bitstr, &size); // lookup
					
					if (size) {			// the corresponding bucket is not empty
						nd += size;
						for (int c = 0; c < size; c++) {
							index = arr[c];
							if (!counter->get(index)) { // if it is not a duplicate
								counter->set(index);
								hammd = match(codes + (UINT64)index*(B_over_8), query, B_over_8);
								nc++;
								if (hammd <= D && numres[hammd] < maxres) {
									res[hammd * K + numres[hammd]] = index + 1;
								}
								numres[hammd]++;
							}
						}
					}
					/* end of processing */

					while (++bit < s && power[bit] == power[bit+1]-1) {
						bitstr ^= (UINT64)1 << (power[bit]-1);
						power[bit] = bit;
					}
					if (bit == s)
						break;
				}
			}

			n = n + numres[s*m+k]; // This line is very tricky ;)
			// The k'th substring (0 based) is the last chance of an
			// item at a Hamming distance of s*m+k to be
			// found. Because if until the k'th substring, an item
			// with distance of s*m+k is not found, then it means that
			// all of the substrings so far have a distance of (s+1)
			// or more, and the remaining substrings have a distance
			// of s or more (total > s*m+k).

			if (n >= maxres)
				break;
		}
	}

	end = clock();

	stats->ticks = end-start;
	stats->numcand = nc;
	stats->numdups = nd;
	stats->numlookups = nl;

	n = 0;
	for (s = 0; s <= D && n < K; s++ ) {
		for (int c = 0; c < numres[s] && n < K; c++)
			results[n++] = res[s*K + c];
	}

	UINT32 total = 0;
	stats->maxrho = -1;
	for (int i=0; i<=B; i++) {
		total += numres[i];
		if (total >= K && stats->maxrho == -1)
			stats->maxrho = i;
	}
	stats->numres = n;
}

void mihasher::cosinequery(UINT32 *results, UINT32* numres, qstat *stats, UINT8 *query, UINT64 * chunks, UINT32 * res, int& max, bool& secondcase,UINT64& num_comp, UINT64& num_buck_check)
{
	UINT32 maxres = K ? K : N;			// if K == 0 that means we want everything to be processed.
	// So maxres = N in that case. Otherwise K limits the results processed.

	UINT32 n = 0; 				// number of results so far obtained (up to a distance of s per chunk)
	UINT32 nc = 0;				// number of candidates tested with full codes (not counting duplicates)
	UINT32 nd = 0;                      	// counting everything retrieved (duplicates are counted multiple times)
	UINT32 nl = 0;				// number of lookups (and xors)
	UINT32 *arr;
	bool checkedtuple[20][20];
	bool flag = false, flag2 = false;
	int size = 0;
	UINT32 numfound = 0;
	UINT32 index;
	int hammd;
	int maxR = 0;
	std::vector<radius> v;
	hammingR HR;
	ofstream resfile,resindex,restuple;
	resfile.open("bucketschecked.txt", ios::out | ios::app);
	resindex.open("NNSindices_64_1M.txt", ios::out | ios::app);
	restuple.open("tuples_64_1M", ios::out | ios::app);
	clock_t start, end;
	radius r,tempr;
	r.r1 = 0;
	r.r2 = 0;
	r.querynorm = 0;
	int norm = 0;
	UINT32 outcome[B/2][B/2][K];
	UINT32 noutcome[B/2][B/2];
	for(int i=0;i<B/2;i++)
		for(int j=0;j<B/2;j++)
		{
			noutcome[i][j] = 0;
			for(int z=0;z<K;z++)
				outcome[i][j][z] = 0;
		}
	for(int i=0;i<B_over_8;i++)
		r.querynorm += popcntll(query[i]);
	norm = r.querynorm;
	int rhat = floor((-1+sqrt(1+4*norm))/2);
	

	start = clock();

	counter->erase();

	split(chunks, query, m, mplus, b);
	
	int s1,s2, r1 ,r2, onesinchunk;
	int curb1k, curb2k;
	UINT64 bitstr2 = 0;
	UINT64 bitstr1 = 0;
	UINT64 bitstr;
	int curb;
	bool coveredR[50];
	int theres[2*K];
	int theresit = 0;

	int bit1,bit2;
	int power1[100], power2[100];
	for(int i=0;i<50;i++)
		coveredR[i] = false;
	coveredR[0] = true;

	for (int i=0;i<20;i++)
		for(int j=0;j<20;j++)
			checkedtuple[i][j] = false;
	UINT64 numchecked=0;
	max = 0;
	while(true)
	{
		if( n>=maxres)
			break;
		int rad = (r.r1+r.r2)/m;
			numchecked += buckets1(B, norm, r.r1, r.r2);
		int in = r.r1>=r.r2 ? 1 : 2; // find the max between r1 and r2
		inttuple candid[rad];
		inttuple e;
		e.x = rad; e.y = 0;
		if(in == 2)
		{
			e.y = rad;
			e.x = 0;
		}
		int cnt = 0;
		for (int z=0;z<rad+1;z++)
		{

			if(e.x<=r.r1 && e.y<=r.r2 && !checkedtuple[e.x][e.y])
			{
				candid[cnt] = e;
				checkedtuple[e.x][e.y] = true;
				cnt++;
			}
			if(in == 1)
			{
				e.x = e.x - 1;
				e.y = e.y + 1;
			}
			else
			{
				e.x = e.x + 1;
				e.y = e.y - 1;
			}
		}


		for(int z=0; z<cnt;z++)
		{ 
			r1 = candid[z].x;r2 = candid[z].y;
			for(int k=0;k<m && numfound<K;k++)
			{
				if (k < mplus)
					curb = b;
				else
					curb = b-1;
				UINT64 chunksk = chunks[k];
				s2= popcntll(chunks[k]);
				if (k < mplus)
					s1 = b-s2;
				else
					s1 = b-1-s2;
				if( r1>s1 || r2>s2)
					continue;
				for (int i=0; i<r1; i++)
					power1[i] = i;			// power[i] stores the location of the i'th 1
				power1[r1] = s1+1;

				bit1 = r1-1;

				int a;
				bitstr1 = 0;
				while(true) // Loop for generating radius of r1
				{

					if(bit1!=-1)
						generatebitstr(power1,bit1,bitstr1);


					else
					{
						for (int i=0; i<r2; i++)
							power2[i] = i;			// power[i] stores the location of the i'th 1
						power2[r2] = s2+1;
						bit2 = r2-1;
						bitstr2 = 0;

						while(true) // Loop for generating radius of r2
						{
							if(bit2!=-1)
								generatebitstr(power2,bit2,bitstr2);

							else
							{

								bitstr = reorderbits(bitstr1,bitstr2,chunksk,curb);

								arr = H[k].query(chunksk ^ bitstr, &size); // lookup
								num_buck_check++;
								if(flag && k==1)
								{
												getchar();
								}
					
								if (size) {			// the corresponding bucket is not empty
										
									nd += size;
									for (int c = 0; c < size; c++) {
										index = arr[c];
										


										if (!counter->get(index)) { // if it is not a duplicate
											counter->set(index);
											num_comp++;
											
											HR = cosinematch(codes + (UINT64)index*(B_over_8), query, B_over_8);
											if(noutcome[HR.r1][HR.r2]<K)
												outcome[HR.r1][HR.r2][noutcome[HR.r1][HR.r2]] = index;
											noutcome[HR.r1][HR.r2]++;
										}
									}
								}
								updatebitstr(power2,bit2,bitstr2,r2);
							}
							if(bit2==r2)
								break;

						}
						updatebitstr(power1,bit1,bitstr1,r1);
						if(bit1==r1)
							break;
					}

				}
			}
		}
		numfound  = numfound + noutcome[r.r1][r.r2];
		if(numfound>=K)
			break;
		tempr.r1 =r.r1;tempr.r2 = r.r2;
		updateR(r,query, coveredR, maxR, v, max);
		if(r.r1+r.r2>rhat)
		  secondcase = true;
	}

}


void mihasher::updateR(radius& R, UINT8* query, bool* coveredR,int& maxR,vector<radius>& v,int& maxheapsize)
{
	int r = R.r1+R.r2;
	if (v.empty() && r==0)
		v.clear();
	if(R.r1==0 && !coveredR[R.r2+1] && R.r2+1<=(B-R.querynorm))
	{
		radius newR;
		newR.r1 = 0;
		newR.r2 = R.r2 + 1;
		newR.querynorm = R.querynorm;
		v.push_back(newR);
		maxR++;
		coveredR[newR.r2] = true;
	}
	if(R.r1<R.querynorm && R.r2!=0)
	{

		radius newR;
		newR.r1 = R.r1 + 1;
		newR.r2 = R.r2 - 1;
		newR.querynorm = R.querynorm;
		v.push_back(newR);
	}
	if(R.r2==0)
	{
		if(!coveredR[r+1] && R.r2+1<=(B-R.querynorm))
		{
			radius newR;
			newR.r1 = 0;
			newR.r2 = r+1;
			newR.querynorm = R.querynorm;
			v.push_back(newR);
			maxR++;
		}
	}
	if(v.empty())
	{
		getchar();
		printf("Vector Empty->Error\n");

	}
	
	if(v.size()>maxheapsize)
	  maxheapsize = v.size();
	std::make_heap(v.begin(),v.end(),comp);
	R = v.front();
	std::pop_heap(v.begin(),v.end(),comp);
	v.pop_back();
}

void mihasher::generatebitstr(int* power, int& bit, UINT64& bitstr)
{
	bitstr ^= (power[bit] == bit) ? (UINT64)1 << power[bit] : (UINT64)3 << (power[bit]-1);
	power[bit]++;
	bit--;


}
void mihasher::updatebitstr(int* power, int& bit, UINT64& bitstr, int s)
{
	while (++bit < s && power[bit] == power[bit+1]-1) {
		bitstr ^= (UINT64)1 << (power[bit]-1);
		power[bit] = bit;
	}
}

mihasher::mihasher(int _B, int _m)
{
	B = _B;
	B_over_8 = B/8;
	m = _m;
	b = ceil((double)B/m);

	D = ceil(B/2.0);		// assuming that B/2 is large enough radius to include all of the k nearest neighbors
	d = ceil((double)D/m);

	mplus = B - m * (b-1);
	// mplus     is the number of chunks with b bits
	// (m-mplus) is the number of chunks with (b-1) bits

	xornum = new UINT32 [d+2];
	xornum[0] = 0;
	for (int i=0; i<=d; i++)
		xornum[i+1] = xornum[i] + choose(b, i);

	H = new SparseHashtable[m];
	// H[i].init might fail
	for (int i=0; i<mplus; i++)
		H[i].init(b);
	for (int i=mplus; i<m; i++)
		H[i].init(b-1);
}




void mihasher::setK(int _K)
{
	K = _K;
}

mihasher::~mihasher()
{
	delete[] xornum;
	delete[] H;
}

void mihasher::binary(UINT64 num)
{
	int rem;

	if (num <= 1)
	{
		cout<<num;
		return;
	}
	rem = num % 2;
	binary(num / 2);
	cout << rem;
}





void mihasher::populate(UINT8 *_codes, UINT32 _N, int dim1codes)
{
	N = _N;
	codes = _codes;

	int k = 0;
	//#pragma omp parallel shared(k)
	{
		UINT64 * chunks = new UINT64[m];
		//#pragma omp for
		for (k=0; k<m; k++) {
			UINT8 * pcodes = codes;
			for (UINT64 i=0; i<N; i++) {
				split(chunks, pcodes, m, mplus, b);
				
				H[k].count_insert(chunks[k], i);
				
				if (i % (int)ceil(N/1000) == 0) {
					printf("%.2f%%\r", (double)i/N * 100);
					fflush(stdout);
				}
				pcodes += dim1codes;
			}


			pcodes = codes;
			for (UINT64 i=0; i<N; i++) {
				split(chunks, pcodes, m, mplus, b);

				H[k].data_insert(chunks[k], i);

				if (i % (int)ceil(N/1000) == 0) {
					printf("%.2f%%\r", (double)i/N * 100);
					fflush(stdout);
				}
				pcodes += dim1codes;
			}
		}

		delete [] chunks;
	}

}
