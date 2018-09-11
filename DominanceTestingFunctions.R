library(primes)
library(Rmpfr)


#Some initial functions required for the dominance testing functions.

#DP
#-Input:  Some integer i in 1-n, and the adjecency martrix of the CP-net structure, A
#-Output: The number of descendent paths in the CP-net, originating at variable i
DP<-function(i,A){
  a=A[i,]
  d=0
  while(sum(a)>0){
    d=d+sum(a)
    a=a%*%A
  }
  d
}

#ancestor
#-Input: Some integer i in 1-n, and the adjecency martrix of the CP-net structure, A
#-Output: An n-length vector which has non-zero jth entry if and only if variable j is an
#         ancestor of variable i
ancestor<-function(i,A){
  a=A[,i]
  b=A[,i]
  while(sum(b)!=0){
    b=b%*%t(A)
    a=a+b
  }
  a
}

#n.val
#-Input: Some integer i in 1-n, and the conditional preference tables (CPTs) of the CP-net
#-Output: The size of the domain of variable i
n.val<-function(i,CPT){
  if(length(dim(CPT[i][[1]]))==0){
    length(CPT[i][[1]])
  }
  else{
    dim(CPT[i][[1]])[length(dim(CPT[i][[1]]))]
  }
}

#Rank
#-Input: Some outcome o, a CP-net (adjacency matrix, A, and CPTs, CPT), level of precision, dig
#-Output: The rank of outcome o, to the specified precision
Rank<-function(o,A,CPT,dig){
  n=dim(A)[1]
  anc=c()
  Anc=list()
  for(i in 1:n){
    a=which(ancestor(i,A)!=0)
    Anc[[i]]=a
    anc=c(anc,length(a))
  }
  RTO=order(-anc)
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  D=mpfrArray(NA,dim = n, precBits = dig)
  for(i in RTO){
    ch=which(A[i,]!=0)
    d=mpfr(0,dig)
    for(j in ch){
      d=mpfr(d+D[j]+1,dig)
    }
    D[i]=d
  }
  rank=mpfr(0,dig)
  for(i in 1:n){
    pa=which(A[,i]!=0)
    pref=CPT[[i]][matrix(c(o[pa],o[i]),1)]
    rank=mpfr(rank+((prod(N[Anc[[i]]]))^(-1)*(D[i]+1)*(N[i]-pref+1)/N[i]),dig)
  }
  return(rank)
}

#Penalty
#-Input: Some outcome o, a CP-net (adjacency matrix, A, and CPTs, CPT), level of precision, dig
#-Output: The penalty of outcome o, to the specified precision
Penalty<-function(o,A,CPT,dig){
  n=dim(A)[1]
  anc=c()
  for(i in 1:n){
    L=length(which(ancestor(i,A)!=0))
    anc=c(anc,L)
  }
  RTO=order(-anc)
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  W=mpfrArray(NA,dim = n, precBits = dig)
  for(i in RTO){
    ch=which(A[i,]!=0)
    w=mpfr(1,dig)
    for(j in ch){
      w=mpfr(w+W[j]*(N[j]-1),dig)
    }
    W[i]=w
  }
  pen=mpfr(0,dig)
  for(i in 1:n){
    pa=which(A[,i]!=0)
    pref=CPT[[i]][matrix(c(o[pa],o[i]),1)]
    d=pref-1
    pen=mpfr(pen+(d*W[i]),dig)
  }
  pen
}

#Pen.Rank
#-Input: Some outcome o, a CP-net (adjacency matrix, A, and CPTs, CPT), level of precision, dig
#-Output: The penalty and rank of outcome o, to the specified precision, as a list
Pen.Rank<-function(o,A,CPT,dig){
  n=dim(A)[1]
  anc=c()
  Anc=list()
  for(i in 1:n){
    a=which(ancestor(i,A)!=0)
    Anc[[i]]=a
    anc=c(anc,length(a))
  }
  RTO=order(-anc)
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  W=mpfrArray(NA,dim = n, precBits = dig)
  D=mpfrArray(NA,dim = n, precBits = dig)
  for(i in RTO){
    ch=which(A[i,]!=0)
    w=mpfr(1,dig)
    d=mpfr(1,dig)
    for(j in ch){
      w=mpfr(w+W[j]*(N[j]-1),dig)
      d=mpfr(d+D[j],dig)
    }
    W[i]=w
    D[i]=d
  }
  pen=mpfr(0,dig)
  rank=mpfr(0,dig)
  for(i in 1:n){
    pa=which(A[,i]!=0)
    pref=CPT[[i]][matrix(c(o[pa],o[i]),1)]
    pen=mpfr(pen+W[i]*(pref-1),dig)
    rank=mpfr(rank+((prod(N[Anc[[i]]]))^(-1)*(D[i])*(N[i]-pref+1)/N[i]),dig)
  }
  return(list(Penalty=pen,Rank=rank))
}


#Dominance Testing Functions

#INPUTS:
#CP-Net consisting of adjacency matrix (A) and CPTs (CPT)
#See CP-net generator script for CP-net input format
#Dominance query "N|= o1>o2?", outcomes o1, o2 are input as vectors (as described in Appendix B of
#Laing et al. (2018) (arxiv.org/abs/1712.08588) and illustrated in the R script for the CP-net generator)


#First Calculate the required precision for the given CP-net
n=dim(A)[1]
N=c()
for(j in 1:n){
  l=length(dim(CPT[[j]]))
  k=dim(CPT[[j]])[l]
  N=c(N,k)
}

dig=max(c((2*ceiling(log(prod(N),base=2))+1),40))

#Dominance query to be answered : N|= o1>o2?
#A - adjacency matrix of the CP-net's, N's, structure
#CPT - conditional preference tables of CP-net, N

#Dominance Testing Function using Rank or Rank + Suffix fixing pruning
#-Inputs: o1,o2,A,CPT as above
#         priority={"rank" - rank prioritisation of leaves,
#                   "rank.diff" - rank + diff prioritisation of leaves}
#         suffix={FALSE - rank pruning, TRUE - rank + suffix fixing pruning}
#         dig - precision, calculated as above

#-Outputs:  list of length 3:
#           1)DQ result, either "False, N does not entail o1 > o2" or "False, N does not entail o1 > o2"
#           2)Outcomes.Considered, defined in Laing et al. (2018) (arxiv.org/abs/1712.08588) as "outcomes traversed"
#           3)Time.Elapsed in seconds until the function terminates (i.e. time taken to answer the query)

DQ.Rank<-function(o1,o2,A,CPT,priority="rank",suffix=TRUE,dig){
  START=as.numeric(Sys.time())
  #START- time the function was started
  if(all(o2==o1)){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Result="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  # If the two outcomes are the same, o1=o2, then the query is trivially false
  n=length(o1)
  #n=number of variables
  r1=Rank(o1,A,CPT,dig=dig)
  r2=Rank(o2,A,CPT,dig=dig)
  #r1,r2 are the ranks of outcomes o1 and o2
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  #N - vector of variable domain sizes
  PA=list()
  CH=list()
  for(i in 1:n){
    PA[[i]]=which(A[,i]!=0)
    CH[[i]]=which(A[i,]!=0)
  }
  #PA - list of the parent sets for each variable
  #CH - list of the children sets for each variable
  Anc=list()
  anc=c()
  for(i in 1:n){
    a=which(ancestor(i,A)!=0)
    Anc[[i]]=a
    anc=c(anc,length(a))
  }
  #Anc - list of the ancestor sets for each variable
  #anc - vector giving the number of ancestors for each variable
  RTO=order(-anc)
  D=mpfrArray(NA,dim = n, precBits = dig)
  for(i in RTO){
    ch=CH[[i]]
    d=mpfr(0,dig)
    for(j in ch){
      d=mpfr(d+D[j]+1,dig)
    }
    D[i]=d
  }
  Least.Inc=mpfrArray(NA,dim = n, precBits = dig)
  for(i in 1:n){
    Inc=mpfr(prod(N[Anc[[i]]])^(-1)*(D[i]+1)/N[i],dig)
    Dec=mpfr(0,dig)
    for(j in CH[[i]]){
      Dec=mpfr(Dec + prod(N[Anc[[j]]])^(-1)*(D[j]+1)*(N[j]-1)/N[j],dig)
    }
    Least.Inc[i]=mpfr(Inc-Dec,dig)
  }
  #Least.Inc - a (Rmpfr) vector giving the Least Rank Improvement for each variable
  DIFF=which(o1!=o2)
  difference=mpfr(sum(Least.Inc[DIFF]),dig)
  if(roundMpfr(r2+difference,precBits=ceiling(dig/2))>roundMpfr(r1,precBits=ceiling(dig/2))){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Result="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  #If r(o1)< r(o2) + L_D(o1,o2) then the query is trivially false
  #(r - rank function, L_D - Least rank difference function)
  if(suffix==TRUE){
    TO=order(anc)
    #TO - a topological ordering of the variables
  }
  count=1
  #count - the number of outcomes in the search tree, currently it is just o2
  if(priority=="rank.diff"){
    S=mpfrArray(o2,precBits=dig)
    S=c(S,r2+difference)
    Search=mpfr2array(S,c(1,(n+1)))
  }
  if(priority=="rank"){
    S=mpfrArray(o2,precBits=dig)
    S=c(S,r2)
    Search=mpfr2array(S,c(1,(n+1)))
  }
  #Search is a matrix with (n+1) columns. This matrix stores the outcomes currently in the search tree.
  #The first n columns of a row give the outcome, o. The end column gives either r(o) or r(o)+L_D(o,o1)
  #depending on whether rank or rank + diff prioritisation is being used.
  #Currently the only outcome in the tree is o2.
  len=1
  #len - the number of outcomes we have considered thus far plus 1
  #Note that we have not yet considered any outcomes
  Q=1
  #Q - the number of unconsidered leaves in the search tree - currently only o2
  while(Q!=0){
    #while there are unconsidered leaves left in the tree:
    o=Search[len,1:n]
    o.v=asNumeric(o)
    #o - next leaf of the tree to be considered (selected by choosing the first listed of the unconsidered leaves in Search)
    #o.v - same outcome as o, with entries stored as numbers
    if(suffix==TRUE){
      ss=max(which(o.v[TO]!=o1[TO]))
      if(ss==n){
        SUFF=integer(0)
      }
      else{
        SUFF=TO[(ss+1):n]
      }
      #SUFF - the variables in the matching suffix of o1 and o
      #(suffix calculated according to topological ordering TO)
    }
    for(i in 1:n){
      if(suffix==TRUE){
        if(!(i%in%SUFF)){
          #If we are using suffix fixing, then for all variables, i, NOT in the matching suffix:
          p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
          #p - preference position of the value taken by variable i in o (given the parents of i take their values in o)
          if(p!=1){
            #If the preference position of variable i can be improved by changing its value:
            pref=c()
            for(j in 1:N[i]){
              pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
            }
            flip=which(pref<p)
            #flip - The values of variable i that would be prefered to the one it takes in o
            for(j in flip){
              #for each of these 'better' values, j:
              o3=o
              o3[i]=j
              o3.v=asNumeric(o3)
              #o3 - outcome obtained from o by flipping the value of variable i to the better value j
              if(all(o3==o1)){
                TIMER=as.numeric(Sys.time())-START
                return(list(DQ.Result="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
              }
              #If o3=o1 then we have reached o1 from o2 by improving flips. Thus the query is true.
              b=FALSE
              if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
                b=TRUE
              }
              if(b){
                #If o3 (the new flip) is not already in the tree:
                r=Rank(o3.v,A,CPT,dig=dig)
                DIFF=which(o1!=o3)
                difference=mpfr(sum(Least.Inc[DIFF]),dig)
                if(roundMpfr(r+difference,precBits=ceiling(dig/2))<=roundMpfr(r1,precBits=ceiling(dig/2))){
                  if(priority=="rank.diff"){
                    k=dim(Search)[1]
                    Search=t(mpfr2array(c(t(Search),o3,(r+difference)),c((n+1),(k+1))))
                  }
                  if(priority=="rank"){
                    k=dim(Search)[1]
                    Search=t(mpfr2array(c(t(Search),o3,r),c((n+1),(k+1))))
                  }
                  #If r(o1)>=r(o3)+L_D(o1,o3), then add o3 to the tree, i.e to Search
                  #(along with the r(o3) or r(o3)+L_D(o1,o3) value depending on prioritisation choice)
                  count=count+1
                  #increase the count of the number of outcomes in the tree by 1 as we have added o3
                }
              }
            }
          }
        }
      }
      else{
        #If not using Suffix fixing, then for each variable, i:
        p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
        #p - preference position of the value taken by variable i in o (given the parents of i take their values in o)
        if(p!=1){
          #If the preference position of variable i can be improved by changing its value:
          pref=c()
          for(j in 1:N[i]){
            pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
          }
          flip=which(pref<p)
          #flip - The values of variable i that would be prefered to the one it takes in o
          for(j in flip){
            #for each of these 'better' values, j:
            o3=o
            o3[i]=j
            o3.v=asNumeric(o3)
            #o3 - outcome obtained from o by flipping the value of variable i to the better value j
            if(all(o3==o1)){
              TIMER=as.numeric(Sys.time())-START
              return(list(DQ.Result="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
            }
            #If o3=o1 then we have reached o1 from o2 by improving flips. Thus the query is true.
            b=FALSE
            if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
              b=TRUE
            }
            if(b){
              #If o3 (the new flip) is not already in the tree:
              r=Rank(o3.v,A,CPT,dig=dig)
              DIFF=which(o1!=o3)
              difference=mpfr(sum(Least.Inc[DIFF]),dig)
              if(roundMpfr(r+difference,precBits=ceiling(dig/2))<=roundMpfr(r1,precBits=ceiling(dig/2))){
                if(priority=="rank.diff"){
                  k=dim(Search)[1]
                  Search=t(mpfr2array(c(t(Search),o3,(r+difference)),c((n+1),(k+1))))
                }
                if(priority=="rank"){
                  k=dim(Search)[1]
                  Search=t(mpfr2array(c(t(Search),o3,r),c((n+1),(k+1))))
                }
                #If r(o1)>=r(o3)+L_D(o1,o3), then add o3 to the tree, i.e to Search
                #(along with the r(o3) or r(o3)+L_D(o1,o3) value depending on prioritisation choice)
                count=count+1
                #increase the count of the number of outcomes in the tree by 1 as we have added o3
              }
            }
          }
        }
      }
    }
    l=dim(Search)[1]
    #l-number of outcomes in the tree currently
    if((l-len)>1){
      #If there are 2+ unconsidered leaves in the tree:
      Search[((len+1):l),]=Search[(order(-Search[((len+1):l),(n+1)])+len),]
      #These leaves will all be at the bottom of Search as we add leaves to the tree by adding them 
      #to the bottom of Search and we pick the next leaf to consider by selecting the top unconsidered leaf
      
      #reorder these unconsidered leaves in Search by orderung them in terms of decreasing 
      #r(*) or r(*)+L_D(*,o1) value (whichever is stored in the n+1 column of Search - recall this is 
      #determined by the prioritisation method chosen)  - This is the implementation of rank/rank+diff prioritisation
    }
    Q=l-len
    #update the number of unconsidered leaves in the tree
    len=len+1
    #o has now been considered so we increase len by 1
  }
  TIMER=as.numeric(Sys.time())-START
  return(list(DQ.Result="False, N does not entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
  #If there are no leaves left to consider (and we did not reach o1) then the search tree is fully constructed
  #and thus the query must be false - we cannot reach o1 from o2 by improving flips
}


#Dominance Testing Function using Penalty or Penalty + Suffix fixing pruning (both with penalty prioritisation of leaves)
#-Inputs: o1,o2,A,CPT as above
#         suffix={FALSE - penalty pruning, TRUE - penalty + suffix fixing pruning}
#         dig - precision, calculated as above

#-Outputs:  list of length 3:
#           1)DQ result, either "False, N does not entail o1 > o2" or "False, N does not entail o1 > o2"
#           2)Outcomes.Considered, defined in Laing et al. (2018) (arxiv.org/abs/1712.08588) as "outcomes traversed"
#           3)Time.Elapsed in seconds until the function terminates (i.e. time taken to answer the query)

DQ.Penalty<-function(o1,o2,A,CPT,suffix=TRUE,dig){
  START=as.numeric(Sys.time())
  #START- time the function was started
  if(all(o2==o1)){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  # If the two outcomes are the same, o1=o2, then the query is trivially false
  p1=Penalty(o1,A,CPT,dig=dig)
  p2=Penalty(o2,A,CPT,dig=dig)
  #p1,p2 are the penalties of outcomes o1 and o2
  f=p2 - length(which(o1!=o2)) - p1
  if(f<0){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  #If the evaluation function, f, is negative, then the query is trivially false
  #(f - the evaluation function based upon penalty values)
  n=length(o1)
  #n=number of variables
  if(suffix==TRUE){
    anc=c()
    for(i in 1:n){
      a=which(ancestor(i,A)!=0)
      anc=c(anc,length(a))
    }
    #anc - vector giving the number of ancestors for each variable
    TO=order(anc)
    #TO - a topological ordering of the variables
  }
  S=mpfrArray(o2,precBits=dig)
  S=c(S,f)
  Search=mpfr2array(S,c(1,(n+1)))
  #Search is a matrix with (n+1) columns. This matrix stores the outcomes currently in the search tree.
  #The first n columns of a row give the outcome, o. The end column gives f(o)
  #Currently the only outcome in the tree is o2.
  count=1
  #count - the number of outcomes in the search tree
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  #N - vector of variable domain sizes
  PA=list()
  for(i in 1:n){
    PA[[i]]=which(A[,i]!=0)
  }
  #PA - list of the parent sets for each variable
  len=1
  #len - the number of outcomes we have considered thus far plus 1
  #Note that we have not yet considered any outcomes
  Q=1
  #Q - the number of unconsidered leaves in the search tree - currently only o2
  while(Q!=0){
    #while there are unconsidered leaves left in the tree:
    o=Search[len,1:n]
    o.v=asNumeric(o)
    #o - next leaf of the tree to be considered (selected by choosing the first listed of the unconsidered leaves in Search)
    #o.v - same outcome as o, with entries stored as numbers
    if(suffix==TRUE){
      ss=max(which(o.v[TO]!=o1[TO]))
      if(ss==n){
        SUFF=integer(0)
      }
      else{
        SUFF=TO[(ss+1):n]
      }
      #SUFF - the variables in the matching suffix of o1 and o
      #(suffix calculated according to topological ordering TO)
    }
    for(i in 1:n){
      if(suffix==TRUE){
        if(!(i%in%SUFF)){
          #If we are using suffix fixing, then for all variables, i, NOT in the matching suffix:
          p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
          #p - preference position of the value taken by variable i in o (given the parents of i take their values in o)
          if(p!=1){
            #If the preference position of variable i can be improved by changing its value:
            pref=c()
            for(j in 1:N[i]){
              pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
            }
            flip=which(pref<p)
            #flip - The values of variable i that would be prefered to the one it takes in o
            for(j in flip){
              #for each of these 'better' values, j:
              o3=o
              o3[i]=j
              o3.v=asNumeric(o3)
              #o3 - outcome obtained from o by flipping the value of variable i to the better value j
              if(all(o3==o1)){
                TIMER=as.numeric(Sys.time())-START
                return(list(DQ.Reult="True, N does entail o1 > o2",Outcomes.Considered=count, Time.Elapsed=TIMER))
              }
              #If o3=o1 then we have reached o1 from o2 by improving flips. Thus the query is true.
              b=FALSE
              if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
                b=TRUE
              }
              if(b){
                #If o3 (the new flip) is not already in the tree:
                f=Penalty(o3.v,A,CPT,dig=dig) - length(which(o1!=o3)) - p1
                if(f>=0){
                  k=dim(Search)[1]
                  Search= t(mpfr2array(c(t(Search),o3,f),c((n+1),(k+1))))
                  #If f(o3)>=0, then add o3 to the tree, i.e to Search, along with the value f(o3)
                  count=count+1
                  #increase the count of the number of outcomes in the tree by 1 as we have added o3
                }
              }
            }
          }
        }
      }
      else{
        #If not using Suffix fixing, then for each variable, i:
        p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
        #p - preference position of the value taken by variable i in o (given the parents of i take their values in o)
        if(p!=1){
          #If the preference position of variable i can be improved by changing its value:
          pref=c()
          for(j in 1:N[i]){
            pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
          }
          flip=which(pref<p)
          #flip - The values of variable i that would be prefered to the one it takes in o
          for(j in flip){
            #for each of these 'better' values, j:
            o3=o
            o3[i]=j
            o3.v=asNumeric(o3)
            #o3 - outcome obtained from o by flipping the value of variable i to the better value j
            if(all(o3==o1)){
              TIMER=as.numeric(Sys.time())-START
              return(list(DQ.Reult="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
            }
            #If o3=o1 then we have reached o1 from o2 by improving flips. Thus the query is true.
            b=FALSE
            if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
              b=TRUE
            }
            if(b){
              #If o3 (the new flip) is not already in the tree:
              f=Penalty(o3.v,A,CPT,dig=dig) - length(which(o1!=o3)) - p1
              if(f>=0){
                k=dim(Search)[1]
                Search= t(mpfr2array(c(t(Search),o3,f),c((n+1),(k+1))))
                #If f(o3)>=0, then add o3 to the tree, i.e to Search, along with the value f(o3)
                count=count+1
                #increase the count of the number of outcomes in the tree by 1 as we have added o3
              }
            }
          }
        }
      }
    }
    l=dim(Search)[1]
    #l-number of outcomes in the tree currently
    if((l-len)>1){
      #If there are 2+ unconsidered leaves in the tree:
      Search[((len+1):l),]=Search[(order(Search[((len+1):l),(n+1)])+len),]
      #These leaves will all be at the bottom of Search as we add leaves to the tree by adding them 
      #to the bottom of Search and we pick the next leaf to consider by selecting the top unconsidered leaf
      
      #reorder these unconsidered leaves in Search by orderung them in terms of increasing f(*) value
      #(stored in the n+1 column of Search) - This is the implementation of penalty prioritisation
    }
    Q=l-len
    #update the number of unconsidered leaves in the tree
    len=len+1
    #o has now been considered so we increase len by 1
  }
  TIMER=as.numeric(Sys.time())-START
  return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
  #If there are no leaves left to consider (and we did not reach o1) then the search tree is fully constructed
  #and thus the query must be false - we cannot reach o1 from o2 by improving flips
}

#Dominance Testing Function using Suffix fixing (with Minimal Depth prioritisation of leaves)
#-Inputs: o1,o2,A,CPT as above
#         dig - precision, calculated as above

#-Outputs:  list of length 3:
#           1)DQ result, either "False, N does not entail o1 > o2" or "False, N does not entail o1 > o2"
#           2)Outcomes.Considered, defined in Laing et al. (2018) (arxiv.org/abs/1712.08588) as "outcomes traversed"
#           3)Time.Elapsed in seconds until the function terminates (i.e. time taken to answer the query)

DQ.SF<-function(o1,o2,A,CPT,dig){
  START=as.numeric(Sys.time())
  #START- time the function was started
  if(all(o2==o1)){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  # If the two outcomes are the same, o1=o2, then the query is trivially false
  n=length(o1)
  #n=number of variables
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  #N - vector of variable domain sizes
  PA=list()
  for(i in 1:n){
    PA[[i]]=which(A[,i]!=0)
  }
  #PA - list of the parent sets for each variable
  anc=c()
  for(i in 1:n){
    a=which(ancestor(i,A)!=0)
    anc=c(anc,length(a))
  }
  #anc - vector giving the number of ancestors for each variable
  TO=order(anc)
  #TO - a topological ordering of the variables
  S=mpfrArray(o2,precBits=dig)
  Search=mpfr2array(S,c(1,(n)))
  #Search is a matrix with (n) columns. This matrix stores the outcomes currently in the search tree.
  #Each row gives an outcome, o, in the tree.
  #Currently the only outcome in the tree is o2.
  count=1
  #count - the number of outcomes in the search tree
  len=1
  #len - the number of outcomes we have considered thus far plus 1
  #Note that we have not yet considered any outcomes
  Q=1
  #Q - the number of unconsidered leaves in the search tree - currently only o2
  while(Q!=0){
    #while there are unconsidered leaves left in the tree:
    o=Search[len,1:n]
    o.v=asNumeric(o)
    #o - next leaf of the tree to be considered (selected by choosing the first listed of the unconsidered leaves in Search)
    #o.v - same outcome as o, with entries stored as numbers
    ss=max(which(o.v[TO]!=o1[TO]))
    if(ss==n){
      SUFF=integer(0)
    }
    else{
      SUFF=TO[(ss+1):n]
    }
    #SUFF - the variables in the matching suffix of o1 and o
    #(suffix calculated according to topological ordering TO)
    for(i in 1:n){
      if(!(i%in%SUFF)){
        #For all variables, i, NOT in the matching suffix:
        p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
        #p - preference position of the value taken by variable i in o (given the parents of i take their values in o)
        if(p!=1){
          #If the preference position of variable i can be improved by changing its value:
          pref=c()
          for(j in 1:N[i]){
            pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
          }
          flip=which(pref<p)
          #flip - The values of variable i that would be prefered to the one it takes in o
          for(j in flip){
            #for each of these 'better' values, j:
            o3=o
            o3[i]=j
            o3.v=asNumeric(o3)
            #o3 - outcome obtained from o by flipping the value of variable i to the better value j
            if(all(o3==o1)){
              TIMER=as.numeric(Sys.time())-START
              return(list(DQ.Reult="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
            }
            #If o3=o1 then we have reached o1 from o2 by improving flips. Thus the query is true.
            b=FALSE
            if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
              b=TRUE
            }
            if(b){
              #If o3 (the new flip) is not already in the tree:
              k=dim(Search)[1]
              Search= t(mpfr2array(c(t(Search),o3),c((n),(k+1))))
              #Add o3 to the tree, i.e to Search
              count=count+1
              #Increase the count of the number of outcomes in the tree by 1 as we have added o3
            }
          }
        }
      }
    }
    l=dim(Search)[1]
    #l-number of outcomes in the tree currently
    Q=l-len
    #update the number of unconsidered leaves in the tree
    len=len+1
    #o has now been considered so we increase len by 1
    
    #Note that here, unlike in the other functions, no re-ordering of the outcomes in Search is required
    #This is because we start with a list of 1 outcome at depth 0. We consider the outcomes in the list 
    #in order and each considered outcome adds any new outcomes to the bottom of the list. All of these 
    #new outcomes will have depth 1 greater than the outcome considered. By this process, the outcome 
    #considered will always have minimal depth out of the unconsidered outcomes (unconsidered leaves).
    #Thus, we are prioritising unconsidered leaves at minimal depth.
  }
  TIMER=as.numeric(Sys.time())-START
  return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
  #If there are no leaves left to consider (and we did not reach o1) then the search tree is fully constructed
  #and thus the query must be false - we cannot reach o1 from o2 by improving flips
}


#Dominance Testing Function using Penalty + Rank or Penalty + Rank + Suffix fixing pruning
#-Inputs: o1,o2,A,CPT as above
#         priority={"rank" - rank prioritisation of leaves,
#                   "rank.diff" - rank + diff prioritisation of leaves,
#                   "penalty" - penalty prioritisation of leaves}
#         suffix={FALSE - Penalty + rank pruning, TRUE - Penalty + rank + suffix fixing pruning}
#         dig - precision, calculated as above

#-Outputs:  list of length 3:
#           1)DQ result, either "False, N does not entail o1 > o2" or "False, N does not entail o1 > o2"
#           2)Outcomes.Considered, defined in Laing et al. (2018) (arxiv.org/abs/1712.08588) as "outcomes traversed"
#           3)Time.Elapsed in seconds until the function terminates (i.e. time taken to answer the query)

DQ.PR<-function(o1,o2,A,CPT,priority="penalty",suffix=TRUE,dig){
  START=as.numeric(Sys.time())
  #START- time the function was started
  if(all(o2==o1)){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  # If the two outcomes are the same, o1=o2, then the query is trivially false
  one=Pen.Rank(o1,A,CPT,dig=dig)
  two=Pen.Rank(o2,A,CPT,dig=dig)
  r1=one$Rank
  r2=two$Rank
  #r1,r2 are the ranks of outcomes o1 and o2
  p1=one$Penalty
  p2=two$Penalty
  #p1,p2 are the penalties of outcomes o1 and o2
  f=p2 - length(which(o1!=o2)) - p1
  if(f<0){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  #If the evaluation function, f, is negative, then the query is trivially false
  #(f - the evaluation function based upon penalty values)
  n=length(o1)
  #n=number of variables
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  #N - vector of variable domain sizes
  PA=list()
  CH=list()
  for(i in 1:n){
    PA[[i]]=which(A[,i]!=0)
    CH[[i]]=which(A[i,]!=0)
  }
  #PA - list of the parent sets for each variable
  #CH - list of the children sets for each variable
  Anc=list()
  anc=c()
  for(i in 1:n){
    a=which(ancestor(i,A)!=0)
    Anc[[i]]=a
    anc=c(anc,length(a))
  }
  #Anc - list of the ancestor sets for each variable
  #anc - vector giving the number of ancestors for each variable
  RTO=order(-anc)
  D=mpfrArray(NA,dim = n, precBits = dig)
  for(i in RTO){
    ch=CH[[i]]
    d=mpfr(0,dig)
    for(j in ch){
      d=mpfr(d+D[j]+1,dig)
    }
    D[i]=d
  }
  Least.Inc=mpfrArray(NA,dim = n, precBits = dig)
  for(i in 1:n){
    Inc=mpfr(prod(N[Anc[[i]]])^(-1)*(D[i]+1)/N[i],dig)
    Dec=mpfr(0,dig)
    for(j in CH[[i]]){
      Dec=mpfr(Dec + prod(N[Anc[[j]]])^(-1)*(D[j]+1)*(N[j]-1)/N[j],dig)
    }
    Least.Inc[i]=mpfr(Inc-Dec,dig)
  }
  #Least.Inc - a vector giving the Least Rank Improvement for each variable
  DIFF=which(o1!=o2)
  difference=mpfr(sum(Least.Inc[DIFF]),dig)
  if( roundMpfr(r2+difference,precBits=ceiling(dig/2))>roundMpfr(r1,precBits=ceiling(dig/2))){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  #If r(o1)< r(o2) + L_D(o1,o2) then the query is trivially false
  #(r - rank function, L_D - Least rank difference function)
  if(suffix==TRUE){
    TO=order(anc)
    #TO - a topological ordering of the variables
  }
  count=1
  #count - the number of outcomes in the search tree, currently it is just o2
  len=1
  #len - the number of outcomes we have considered thus far plus 1
  #Note that we have not yet considered any outcomes
  if(priority=="rank.diff"){
    S=mpfrArray(o2,precBits=dig)
    S=c(S,r2+difference)
    Search=mpfr2array(S,c(1,(n+1)))
  }
  if(priority=="rank"){
    S=mpfrArray(o2,precBits=dig)
    S=c(S,r2)
    Search=mpfr2array(S,c(1,(n+1)))
  }
  if(priority=="penalty"){
    S=mpfrArray(o2,precBits=dig)
    S=c(S,f)
    Search=mpfr2array(S,c(1,(n+1)))
  }
  #Search is a matrix with (n+1) columns. This matrix stores the outcomes currently in the search tree.
  #The first n columns of a row give the outcome, o. The end column gives either f(o), r(o) or 
  #r(o)+L_D(o,o1) depending on whether penalty, rank or rank + diff prioritisation is being used.
  #Currently the only outcome in the tree is o2.
  Q=1
  #Q - the number of unconsidered leaves in the search tree - currently only o2
  while(Q!=0){
    #while there are unconsidered leaves left in the tree:
    o=Search[len,1:n]
    o.v=asNumeric(o)
    #o - next leaf of the tree to be considered (selected by choosing the first listed of the unconsidered leaves in Search)
    #o.v - same outcome as o, with entries stored as numbers
    if(suffix==TRUE){
      ss=max(which(o.v[TO]!=o1[TO]))
      if(ss==n){
        SUFF=integer(0)
      }
      else{
        SUFF=TO[(ss+1):n]
      }
      #SUFF - the variables in the matching suffix of o1 and o
      #(suffix calculated according to topological ordering TO)
    }
    for(i in 1:n){
      if(suffix==TRUE){
        if(!(i%in%SUFF)){
          #If we are using suffix fixing, then for all variables, i, NOT in the matching suffix:
          p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
          #p - preference position of the value taken by variable i in o (given the parents of i take their values in o)
          if(p!=1){
            #If the preference position of variable i can be improved by changing its value:
            pref=c()
            for(j in 1:N[i]){
              pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
            }
            flip=which(pref<p)
            #flip - The values of variable i that would be prefered to the one it takes in o
            for(j in flip){
              #for each of these 'better' values, j:
              o3=o
              o3[i]=j
              o3.v=asNumeric(o3)
              #o3 - outcome obtained from o by flipping the value of variable i to the better value j
              if(all(o1==o3)){
                TIMER=as.numeric(Sys.time())-START
                return(list(DQ.Reult="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
              }
              #If o3=o1 then we have reached o1 from o2 by improving flips. Thus the query is true.
              b=FALSE
              if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
                b=TRUE
              }
              if(b){
                #If o3 (the new flip) is not already in the tree:
                imp=Pen.Rank(o3.v,A,CPT,dig=dig)
                r=imp$Rank
                p=imp$Penalty
                f=p - length(which(o1!=o3)) - p1
                #r is the rank of o3, r(o3)
                #p is the penalty of o3
                #f is the evaluation function value at o3, f(o3)
                if(f>=0){
                  #If f(o3)>=0 :
                  DIFF=which(o1!=o3)
                  difference=mpfr(sum(Least.Inc[DIFF]),dig)
                  if(roundMpfr(r+difference,precBits=ceiling(dig/2))<=roundMpfr(r1,precBits=ceiling(dig/2))){
                    #If r(o1)>=r(o3)+L_D(o1,o3) :
                    if(priority=="rank.diff"){
                      k=dim(Search)[1]
                      Search=t(mpfr2array(c(t(Search),o3,(r+difference)),c((n+1),(k+1))))
                    }
                    if(priority=="rank"){
                      k=dim(Search)[1]
                      Search=t(mpfr2array(c(t(Search),o3,r),c((n+1),(k+1))))
                    }
                    if(priority=="penalty"){
                      k=dim(Search)[1]
                      Search=t(mpfr2array(c(t(Search),o3,f),c((n+1),(k+1))))
                    }
                    #Add o3 to the tree, i.e to Search (along with the f(o3), r(o3) or r(o3)+L_D(o1,o3)
                    #value depending on prioritisation choice)
                    count=count+1
                    #increase the count of the number of outcomes in the tree by 1 as we have added o3
                  }
                }
              }
            }
          }
        }
      }
      else{
        #If not using Suffix fixing, then for each variable, i:
        p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
        #p - preference position of the value taken by variable i in o (given the parents of i take their values in o)
        if(p!=1){
          #If the preference position of variable i can be improved by changing its value:
          pref=c()
          for(j in 1:N[i]){
            pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
          }
          flip=which(pref<p)
          #flip - The values of variable i that would be prefered to the one it takes in o
          for(j in flip){
            #for each of these 'better' values, j:
            o3=o
            o3[i]=j
            o3.v=asNumeric(o3)
            #o3 - outcome obtained from o by flipping the value of variable i to the better value j
            if(all(o1==o3)){
              TIMER=as.numeric(Sys.time())-START
              return(list(DQ.Reult="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
            }
            #If o3=o1 then we have reached o1 from o2 by improving flips. Thus the query is true.
            b=FALSE
            if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
              b=TRUE
            }
            if(b){
              #If o3 (the new flip) is not already in the tree:
              imp=Pen.Rank(o3.v,A,CPT,dig=dig)
              r=imp$Rank
              p=imp$Penalty
              f=p - length(which(o1!=o3)) - p1
              #r is the rank of o3, r(o3)
              #p is the penalty of o3
              #f is the evaluation function value at o3, f(o3)
              if(f>=0){
                #If f(o3)>=0 :
                DIFF=which(o1!=o3)
                difference=mpfr(sum(Least.Inc[DIFF]),dig)
                if(roundMpfr(r+difference,precBits=ceiling(dig/2))<=roundMpfr(r1,precBits=ceiling(dig/2))){
                  #If r(o1)>=r(o3)+L_D(o1,o3) :
                  if(priority=="rank.diff"){
                    k=dim(Search)[1]
                    Search=t(mpfr2array(c(t(Search),o3,(r+difference)),c((n+1),(k+1))))
                  }
                  if(priority=="rank"){
                    k=dim(Search)[1]
                    Search=t(mpfr2array(c(t(Search),o3,r),c((n+1),(k+1))))
                  }
                  if(priority=="penalty"){
                    k=dim(Search)[1]
                    Search=t(mpfr2array(c(t(Search),o3,f),c((n+1),(k+1))))
                  }
                  #Add o3 to the tree, i.e to Search (along with the f(o3), r(o3) or r(o3)+L_D(o1,o3)
                  #value depending on prioritisation choice)
                  count=count+1
                  #increase the count of the number of outcomes in the tree by 1 as we have added o3
                }
              }
            }
          }
        }
      }
    }
    l=dim(Search)[1]
    #l-number of outcomes in the tree currently
    if((l-len)>1){
      #If there are 2+ unconsidered leaves in the tree:
      if(priority=="penalty"){
        Search[((len+1):l),]=Search[(order(Search[((len+1):l),(n+1)])+len),]
      }
      else{
        Search[((len+1):l),]=Search[(order(-Search[((len+1):l),(n+1)])+len),]
      }
      #These leaves will all be at the bottom of Search as we add leaves to the tree by adding them 
      #to the bottom of Search and we pick the next leaf to consider by selecting the top unconsidered leaf
      
      #If we are using penalty prioritisation, reorder these unconsidered leaves in Search by ordering them 
      #in terms of increasing f(*) values (recall, this is stored in the n+1 column of Search)
      #If we are using rank/rank+diff prioritisation, reorder these unconsidered leaves in Search by ordering
      #them in terms of decreasing r(*) or r(*)+L_D(*,o1) value respectively (This is stored in the n+1 column 
      #of Search - recall this is determined by the prioritisation method chosen)
      #- This is the implementation of the prioritisation method chosen
    }
    Q=l-len
    #update the number of unconsidered leaves in the tree
    len=len+1
    #o has now been considered so we increase len by 1
  }
  TIMER=as.numeric(Sys.time())-START
  return(list(DQ.Reult="False, N does not entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
  #If there are no leaves left to consider (and we did not reach o1) then the search tree is fully constructed
  #and thus the query must be false - we cannot reach o1 from o2 by improving flips
}


