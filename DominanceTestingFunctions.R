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
#Dominance query "N|= o1>o2?", outcomes o1, o2 are input as vectors (as described in Appendix B of ARXIV LINK)


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

#Dominance Testing Function for Rank or Rank + Suffix fixing pruning
#-Inputs: o1,o2,A,CPT as above
#         priority={"r" - rank prioritisation, "r.diff" - rank + diff prioritisation}
#         suffix={FALSE - rank pruning, TRUE - rank + suffix fixing pruning}
#         dig - precision, calculated as above

#-Outputs:  list of length 3:
#           1)DQ result, either "False, N does not entail o1 > o2" or "False, N does not entail o1 > o2"
#           2)Outcomes.Considered, as defined in ARXIV LINK
#           3)Time.Elapsed in seconds until the function terminates (i.e. time taken to answer the query)

TE_DQ.Rank.2<-function(o1,o2,A,CPT,priority="r",suffix=TRUE,dig){
  START=as.numeric(Sys.time())
  if(all(o2==o1)){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Result="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  # If the two functions are the same, o1=o2, then the query is trivially false
  n=length(o1)
  r1=Rank(o1,A,CPT,dig=dig)
  r2=Rank(o2,A,CPT,dig=dig)
  N=c()
  for(i in 1:n){
    N=c(N,n.val(i,CPT))
  }
  PA=list()
  CH=list()
  for(i in 1:n){
    PA[[i]]=which(A[,i]!=0)
    CH[[i]]=which(A[i,]!=0)
  }
  ## Finding the vector of least poss increases 
  Anc=list()
  anc=c()
  for(i in 1:n){
    a=which(ancestor(i,A)!=0)
    Anc[[i]]=a
    anc=c(anc,length(a))
  }
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
  ##
  DIFF=which(o1!=o2)
  difference=mpfr(sum(Least.Inc[DIFF]),dig)
  if(roundMpfr(r2+difference,precBits=ceiling(dig/2))>roundMpfr(r1,precBits=ceiling(dig/2))){
    TIMER=as.numeric(Sys.time())-START
    return(list(DQ.Result="False, N does not entail o1 > o2",Outcomes.Considered=0,Time.Elapsed=TIMER))
  }
  if(suffix==TRUE){
    TO=order(anc)
  }
  count=1
  if(priority=="r.diff"){
    S=mpfrArray(o2,precBits=dig)
    S=c(S,r2+difference)
    Search=mpfr2array(S,c(1,(n+1)))
  }
  if(priority=="r"){
    S=mpfrArray(o2,precBits=dig)
    S=c(S,r2)
    Search=mpfr2array(S,c(1,(n+1)))
  }
  len=1
  Q=1
  while(Q!=0){
    o=Search[len,1:n]
    o.v=asNumeric(o)
    if(suffix==TRUE){
      ss=max(which(o.v[TO]!=o1[TO]))
      if(ss==n){
        SUFF=integer(0)
      }
      else{
        SUFF=TO[(ss+1):n]
      }
    }
    for(i in 1:n){
      if(as.numeric(Sys.time())> STOPTIME){
        TIMER=as.numeric(Sys.time())-START
        i<-as.numeric(args[1])
        #i number 1-100 of CP_nets
        j<-as.numeric(args[2])
        #j number 1-10 of DQ
        n<-as.numeric(args[3])
        ccc=n-1
        d<-as.numeric(args[4])
        fun="RANK_RDIFF"
        DATA = list(i,j, fun, o, len, count, TIMER, Search)
        save(DATA,file=paste("/nobackup/mm12kl/ARC_Transfer/data/TE_CPN/TE_INCOMPLETE/MID_FUN_TE_RANK_RDIFF.",n,ccc,d,".",i,j,".RData",sep=""))
        Rank.Diff.r.diff.F.Count=c(Rank.Diff.r.diff.F.Count,NA)
        Rank.Diff.r.diff.F.Time=c(Rank.Diff.r.diff.F.Time,NA)
        Result=c(Result,NA)
        
        save(Result, Rank.Diff.r.diff.F.Count, Rank.Diff.r.diff.F.Time, file=paste("/nobackup/mm12kl/ARC_Transfer/data/TE_CPN/Rank.Diff.r.diff.F.",n,ccc,d,".",i,".RData",sep=""))
        Boolean = FALSE
        stop("BREAK")
      }
      if(suffix==TRUE){
        if(!(i%in%SUFF)){
          p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
          if(p!=1){
            pref=c()
            for(j in 1:N[i]){
              pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
            }
            flip=which(pref<p)
            for(j in flip){
              o3=o
              o3[i]=j
              o3.v=asNumeric(o3)
              if(all(o3==o1)){
                TIMER=as.numeric(Sys.time())-START
                return(list(DQ.Result="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
              }
              b=FALSE
              if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
                b=TRUE
              }
              if(b){
                r=Rank(o3.v,A,CPT,dig=dig)
                DIFF=which(o1!=o3)
                difference=mpfr(sum(Least.Inc[DIFF]),dig)
                if(roundMpfr(r+difference,precBits=ceiling(dig/2))<=roundMpfr(r1,precBits=ceiling(dig/2))){
                  if(priority=="r.diff"){
                    k=dim(Search)[1]
                    Search=t(mpfr2array(c(t(Search),o3,(r+difference)),c((n+1),(k+1))))
                  }
                  if(priority=="r"){
                    k=dim(Search)[1]
                    Search=t(mpfr2array(c(t(Search),o3,r),c((n+1),(k+1))))
                  }
                  count=count+1
                }
              }
            }
          }
        }
      }
      else{
        p=CPT[[i]][matrix(c(o.v[PA[[i]]],o.v[i]),1)]
        if(p!=1){
          pref=c()
          for(j in 1:N[i]){
            pref=c(pref,CPT[[i]][matrix(c(o.v[PA[[i]]],j),1)])
          }
          flip=which(pref<p)
          for(j in flip){
            o3=o
            o3[i]=j
            o3.v=asNumeric(o3)
            if(all(o3==o1)){
              TIMER=as.numeric(Sys.time())-START
              return(list(DQ.Result="True, N does entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
            }
            b=FALSE
            if(!(any(unlist(apply(Search,1,function(x) all.equal(x[1:n],o3)==TRUE))))){
              b=TRUE
            }
            if(b){
              r=Rank(o3.v,A,CPT,dig=dig)
              DIFF=which(o1!=o3)
              difference=mpfr(sum(Least.Inc[DIFF]),dig)
              if(roundMpfr(r+difference,precBits=ceiling(dig/2))<=roundMpfr(r1,precBits=ceiling(dig/2))){
                if(priority=="r.diff"){
                  k=dim(Search)[1]
                  Search=t(mpfr2array(c(t(Search),o3,(r+difference)),c((n+1),(k+1))))
                }
                if(priority=="r"){
                  k=dim(Search)[1]
                  Search=t(mpfr2array(c(t(Search),o3,r),c((n+1),(k+1))))
                }
                count=count+1
              }
            }
          }
        }
      }
    }
    l=dim(Search)[1]
    if((l-len)>1){
      Search[((len+1):l),]=Search[(order(-Search[((len+1):l),(n+1)])+len),]
    }
    Q=l-len
    len=len+1
  }
  TIMER=as.numeric(Sys.time())-START
  return(list(DQ.Result="False, N does not entail o1 > o2",Outcomes.Considered=count,Time.Elapsed=TIMER))
}
