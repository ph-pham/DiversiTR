#include <Rcpp.h>
#include <ctime> 
#include <cmath>        // std::abs


using namespace Rcpp;


// [[Rcpp::export]]
int indiceRand(Rcpp::NumericVector row, int  sumrow){
	int rnd = rand()%sumrow;
	for(int i=0; i<row.size(); i++) {
		if(rnd < row(i))
			return i;
		rnd -= row(i);
	}
}

Rcpp::NumericVector cumsum1(Rcpp::NumericVector x){
  // initialize an accumulator variable
  int acc = 0;
  
  // initialize the result vector
  Rcpp::NumericVector tmp(x.size());
  Rcpp::NumericVector res;

  Rcpp::NumericVector::iterator it = tmp.begin();

  int i =0,j = 0;
  while(i<x.size()){
    if(x[i]>1){
      // std::cout << "plus1 : ";
      // std::cout << x[i];
      // std::cout << "\n";
    }
    if(x[i]>0){
      acc += x[i];
      tmp[j] = acc;
      it++;
      j++;
    }
    i++;
  }
  res.assign(tmp.begin(), it-1);
  // for(int i = 0; i < x.size(); i++){
  //   acc += x[i];
  //   res[i] = acc;
  // }
  return res;
}


float mean(Rcpp::NumericVector lspecies, int sample){;
  int out = 0;
  for(int i = 0; i<sample; i++){
    out+=lspecies(i);
  }
  // std::cout << "out : ";
  // std::cout << out;
  // std::cout << "\n";
  // std::cout << "sample : ";
  // std::cout << sample;
  // std::cout << "\n";
  
  return float(out)/float(sample);
}

Rcpp::NumericVector clear(Rcpp::NumericVector s){
  int len = s.size();
  for(int i = 0; i<len; i++){
    if(s(i)==0)
      return s;
    s(i)=0;
  }
  return s;
}

int dichotomie(Rcpp::NumericVector row, int  rand){
  int size = row.size();
  int i;
  int left = 0;
  int right = rand-1;
  if(rand>=size)
    right = size -1;
  // std::cout << "rand : ";
  // std::cout << rand;
  // std::cout << "\n";
  
  while(1){ 
    i = (left+right)/2;
    Rcpp::checkUserInterrupt();
    
    if(row(i)<=rand){
      if(i+1==size || row(i+1)>=rand){
        // std::cout << "i : ";
        // std::cout << rand;
        // std::cout << "\t";
        // std::cout << row(i);
        // std::cout << "\n";
        return i+1;
      }

      left = i;
      }
    else{
      if(i == 0 || row(i-1)<=rand){
        // std::cout << "i : ";
        // std::cout << rand;
        // std::cout << "\t";
        // std::cout << row(i);
        // std::cout << "\n";
        return i;
      }
      right = i;
      }
    }
}


int indiceRandcumul(Rcpp::NumericVector row, int  sumrow){
  int rnd = rand()%sumrow;
  return dichotomie(row, rnd);
}

//' compute rarefaction
//'
//' This function returns
//' @param df a data frame of counts with counts for each sample in columns and clonotypes in rows
//' @return a dataframe 
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame rarefaction (Rcpp::DataFrame df, int sample){//, int step){
  Rcpp::NumericMatrix x = internal::convert_using_rfunction(df, "as.matrix");

  srand((int)time(0));
  int ncols = x.ncol(), outind, samplesize, sumrow;
  Rcpp::NumericVector sumrows(ncols);
  Rcpp::NumericVector v, vpos, vcumsum;
  Rcpp::NumericVector lspecies(sample);
  for(int i = 0; i<ncols; i++){
	v = x( _, i);
	sumrows(i) = sum(v);//sum(v1);
  }
  int outcols = max(sumrows);
  int step = min(sumrows)/20;
  Rcpp::NumericVector s(outcols);
  std::cout << "step : ";
  std::cout << step;
  std::cout << "\n";
  //cpp::StringVector clonotypes = x(_, 0);
  Rcpp::NumericMatrix out(ceil(float(outcols)/float(step)), ncols+1);
  for (int i = 0; i < ncols; i++){
    //std::cout << "i : ";
    //std::cout << i;
    //std::cout << "\n";
  
	v = x(_, i);
  //vpos = v[v>0];
 
  vcumsum = cumsum1(v);
  std::cout << "size cumsum\n";
  std::cout << vcumsum.size();
  std::cout << "\n";
	sumrow = sumrows(i);
	out(0, i) = 1;
	out(0, ncols)=1;
	outind = 0;
	samplesize = 0;
	//int sumcol = std::accumulate(v1.begin(), v1.end(), 0.0);
	//for( int samplesize = step; samplesize < sumrow; samplesize+=step){
	do{
	  outind++;
	  if(samplesize +step < sumrow)
	    samplesize+=step;
	  else
	    samplesize = sumrow;
	  
	  if(sumrow==outcols){
	    out(outind, ncols) = samplesize;
	  }
	  //std::cout << "samplesize : ";
	  //std::cout << samplesize ;
	  //std::cout << "\n";
	  Rcpp::checkUserInterrupt();
		for (int j = 0; j < sample ; j++){
			for(int k = 0; k < samplesize; k++){
			  //std::cout << "outcol/100 : ";
			  //std::cout << outcols/100 ;
			  //std::cout << "\n";
			  //std::cout << "k : ";
			  //std::cout << k ;
			  //std::cout << "in\n";
			  int ind = indiceRandcumul(vcumsum, sumrow);
				s(k) = ind;

			}
			
			lspecies(j) = unique(s).size();
		  if(samplesize<outcols-1)
		    lspecies(j)--;
		  // std::cout << "lsppecies : ";
		  // std::cout << lspecies ;
		  // std::cout << "\n";
		}
		//std::cout << "coucou\n";
		//free(s);
		out(outind, i) = mean(lspecies, sample);
		std::cout << out(outind, i) ;
		std::cout << "\n";
		
		std::cout << std::abs(out(outind, i) - out(outind-1, i));
		std::cout << "\n";
		
	}while(samplesize <= sumrow && std::abs(out(outind, i) - out(outind-1, i)) >= 1);
	std::fill(s.begin(), s.end(), 0);
	
  }
  Rcpp::StringVector cols = colnames(x);
  cols.push_back("size");
  colnames(out) = cols;
  return wrap(out);
}
