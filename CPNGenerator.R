library(primes)
library(Rmpfr)

#rand.dag
#Inputs - n, the number of variables in the CP-net, p, the max number of parents for any variable
#Output - L, a (list) valid dagcode that corresponds to a DAG over n nodes, each of which has at most
#         p parents

rand.dag<-function(n,p){
  L=list()
  #L is an empty list that we will build up to be the dagcode
  var=1:n
  #var the set of variables in the CP-net, enumerated 1-n
  for(i in 1:(n-1)){
    U=c()
    if(i>1){
      for(j in 1:(i-1)){
        U=union(U,L[[j]])
      }
    }
    #First set U to be the union of all sets of variables already in L
    Bool=F
    while(Bool==F){
      s=sample(0:p,1)
      #s -  a random number between 0 and p
      t=sample(var,s)
      #t - random set of s variables
      if(length(union(U,t))<=i){
        Bool=T
      }
      #If there are no more than i variables in the union of t with all the other sets currently 
      #in L (the dagcode), exit the while loop, othewise return to the start of the while loop
    }
    L[[i]]=t
    #set t to be the next entry of the list, L
  }
  return(L)
  #Output the completed dagcode
}

#make.dag
#Inputs -L, a (list) valid dagcode (recall, this is a list of sets)
#Outputs -A, the adjacency matrix of the DAG that corresponds to dagcode L

make.dag<-function(L){
  n=length(L)+1
  #n - number of variables in the CP-net (i.e. nodes in the DAG)
  A=matrix(rep(0,n^2),ncol=n)
  #A - n by n matrix of 0s (will be the adjacency matrix)
  var=1:n
  #var the set of unassigned variables/nodes
  #currently this is all of the variables (enumerated 1-n)
  for(j in (n-1):1){
    #Working backwards through the entries of the dagcode L:
    #Current entry being considered is L[[j]]
    U=c()
    for(k in 1:j){
      U=union(U,L[[k]])
    }
    #U - the union of L[[j]] with all the sets that come before it in the dagcode
    R=setdiff(var,U)
    i=max(R)
    #i - the largest unassigned variable number that is not in this union
    for(l in L[[j]]){
      A[l,i]=1
    }
    #Set the variables in the current dagcode entry, L[[j]], to be the parents of i
    #Add these parent-child edged to the adjacency matrix
    w=which(var==i)
    var=var[-w]
    #remove i from the set of unnassigned variables as we have assigned parents to i
  }
  return(A)
  #Output the now complete DAG adjancency matrix
}

#encode
#Input - x, Any finite vector of positive integers, dig, the required level of precision.
#Output - An integer that encodes x, i.e. given this output one can recover the initial vector
#         (thus, any two distrinct inputs given distinct encodings)
#Note: This encoding is stored and returned as an Rmpfr object, as this allows us to specify the precision
#      to which it is stored. Otherwise, rounding may occur which means the encodings are no longer unique

encode<-function(x,dig){
  n=length(x)
  #n - length of input vector x
  E=mpfr(1,dig)
  #E=1
  counter=1
  Pr=c()
  while(length(Pr)<n){
    Next=generate_primes(min=10^(counter-1),max=10^counter)
    Pr=c(Pr,Next)
    counter=counter+1
  }
  Pr=Pr[1:n]
  #Pr - vector of the first n primes in order
  for(i in 1:n){
    #For each entry of x:
    #current entry of x is x[i]
    E=mpfr(E*Pr[i]^(x[i]),dig)
    #Multiply E by the corresponding prime entry of Pr, p, to the power x[i]
  }
  return(E)
  #Output E
}

#Deg
#Input - CPT (stored as a list), a conditional prefernce table that may or may not be valid
#            (let us say it is the CPT of variable X)
#        dig, the required level of precision
#Output - {"degenerate" or "non-degenerate"} - Whether the CPT is degenerate or not (not valid or valid)

Deg<-function(CPT,dig){
  p=length(dim(CPT))-1
  #p - the number of parents of X (X-the variable this CPT belongs to)
  d=dim(CPT)[p+1]
  #d - the size of X's domain
  if(p==0){
    return("non-degenerate")
  }
  #If X has no parents, it is not possible for the CPT to be degenerate
  if(p==1){
    #If X has one parent:
    x=c()
    for(j in 1:d){
      o=c(1,j)
      x=c(x,CPT[matrix(o,1)])
    }
    #x - the preference order over X when its parent is set to value 1
    #Note - all variable values are input simply as the enumeration 1 to size of domain
    E=encode(x,dig)
    #E- integer encoding of the preference order x
    for(i in 2:(dim(CPT)[1])){
      #for all other (not 1) possible values the parent can take, i:
      x=c()
      for(j in 1:d){
        o=c(i,j)
        x=c(x,CPT[matrix(o,1)])
      }
      #x - the preference order over X when its parent is set to value i
      T=encode(x,dig)
      #T- integer encoding of the preference order x
      if(T!=E){
        return("non-degenerate")
      }
      #If T=/=E then X has a different preference order under parent=1 and parent=j
      #thus, the only parent is valid and so the CPT is non-degenerate
    }
    return("degenerate")
    #If all of the other parent values give the same preference order as parent=1, then changing the
    #parent's value does not affect preference over X. Thus, the parent is invalid. 
    #So the CPT is degenerate
  }
  #If there are >1 parents of X:
  for(i in 1:p){
    #for each parent of X, i:
    q=c()
    for(j in 1:p){
      if(j!=i){
        q=c(q,(dim(CPT)[j]))
      }
    }
    #q - list of domain sizes for each parent that is not i
    asst=matrix(rep(1,prod(q)*(p-1)),ncol=p-1)
    for(j in 1:(p-1)){
      #d=q[j]
      if(j==(p-1)){
        a=1
      }
      else{
        a=prod(q[(j+1):(p-1)])
      }
      b=prod(q[j:(p-1)])
      b=dim(asst)[1]/b
      asst[,j]=rep(1:q[j],each=a,b)
    }
    #asst is a matrix with (p-1) columns. Each row of asst gives a possible assignment of values to
    #the of parents of X other than i. Each possible assignment of values to the not-i parents is
    #present as exactly one row of asst
    for(j in 1:dim(asst)[1]){
      #For each row of asst, i.e. for each possible assignment of values to the other (not i) parents:
      #call the assignment being considered currently, assignment j
      if(i==1){
        Pa=c(1,asst[j,])
      }
      if(i==p){
        Pa=c(asst[j,],1)
      }
      if(i>1&&i<p){
        Pa=c(asst[j,1:(i-1)],1,asst[j,i:(p-1)])
      }
      #Pa - assignment of values to the parents of X that sets parent i to equal 1, all other parents
      #are assigned the value they take in assignment j
      x=c()
      for(l in 1:d){
        O=c(Pa,l)
        x=c(x, CPT[matrix(O,1)])
      }
      #x - the preference order over X given that parents=Pa
      E=encode(x,dig)
      #E - integer encoding of x
      Bool=0
      for(k in 2:dim(CPT)[i]){
        #For all other possible values parent i could take (other than 1), k:
        if(i==1){
          Pa=c(k,asst[j,])
        }
        if(i==p){
          Pa=c(asst[j,],k)
        }
        if(i>1&&i<p){
          Pa=c(asst[j,1:(i-1)],k,asst[j,i:(p-1)])
        }
        #Pa - assignment of values to the parents of X that sets parent i to equal k, all other parents
        #are assigned the value they take in assignment j
        x=c()
        for(l in 1:d){
          O=c(Pa,l)
          x=c(x, CPT[matrix(O,1)])
        }
        #x - the preference order over X given that parents=Pa
        T=encode(x,dig)
        #T - integer encoding of x
        if(T!=E){
          Bool=1
          break
          #If T=/=E then X has a different preference order under {parent i=1 and other parents=assignment j}
          #and {parent i=k and other parents=assignment j}. Thus, changing only parent i can change
          #preference over X. This confirms that parent i is valid - we do not need to look at the other
          #values we can assign parent i, or the other possible assignments to the other parents.
          #In this case, set bool=1 and break out of the loop considering the other values for parent i
        }
        
        #If preference over X is again the same for {parent i=k and other parents=assignment j} (T=E),
        #then we have not yet shown parent i to be valid.
      }
      if(Bool==1){
        break
        #If we have shown parent i to be valid, we do not need to consider the other possible assignments
        #to the other parents, we can move on to assessing whether the next parent is valid. Thus, we break
        #out of this loop (considering each possible assignment to the other pents)
      }
      #If Bool=0, we have not yet shown parent i to be valid (i.e that changing only parent i can change
      #preference over X. Thus, we move on to considering what happens when we change the value of parent i
      #with the next assignment of values to the other parents
    }
    #We have now either shown parent i to be valid (bool=1) or we have found that no pair of parent assignments
    #that differ only on parent i have different associated preferences over X, thus Bool remains =0
    if(Bool==0){
      return("degenerate")
      #If Bool=0 then parent i is not a valid parent, so the CPT is degenerate
    }
    #If parent i is valid, we must go on to check the next parent - all parents must be valid for the CPT
    #to be non-degenerate
  }
  return("non-degenerate")
  #If all parents have been found valid, then the CPT is non-degenerate
}

#rand.CPT
#Input - A, adjacency matrix of CP-net structure
#        d, maximum domain size
#        dig, required level of precision
#Output - A random list of valid (non-degenerate) CPTs for this CP-net

rand.CPT<-function (A,d,dig){
  n=dim(A)[1]
  #n-total number of variables
  if(d>2){
    N=sample(2:d,n,replace = TRUE)
  }
  else{
    N=rep(2,n)
  }
  #N - vector of domain sizes for the variables, randomly generated
  #Note, each variable must have domain size between 2-d
  CPT=list()
  #CPT - empty list, will be the output list of CPTs
  PA=list()
  for(i in 1:n){
    W=which(A[,i]!=0)
    PA[[i]]=W
  }
  #PA - list of the parent sets for each variable
  DEGENERATE=TRUE
  while(DEGENERATE==TRUE){
    for(i in 1:n){
      #for each CP-net variable, i:
      if(length(PA[[i]])==0){
        CPT[[i]]=array(dim=N[i])
        #If variable i has no parents, CPT(i) is a |Dom(i)|-length vector
        CPT[[i]][1:(N[i])]=sample(N[i])
        #The entries of this CPT is simply a preference order over the values of variable i.
        #The way this is input should be that the kth entry=preference position of kth value of variable i
        #So if the kth entry is 1 then variable i=k is the best choice for variable i. If the entry is 2
        #then its the second best choice. 
        #We generate this preference order by randomly generating an ordering the numbers 1-|Dom(i)|.
      }
      else{
        #If variable i has parents:
        N.Pa=N[PA[[i]]]
        #N.Pa - vector of parent domain sizes
        N.i=c(N.Pa,N[i])
        #N.i - N.Pa with |Dom(i)| at the end
        CPT[[i]]=array(dim = N.i)
        #N.i is the dimension of the CPT(i) array
        o=rep(1,length(PA[[i]]))
        pref=sample(N[i])
        for(j in 1:N[i]){
          o.f=c(o,j)
          CPT[[i]][matrix(o.f,1)]=pref[j]
        }
        #Input to CPT(i) a random ordering of the numbers 1-|Dom(i)| as the preference ordering of i
        #Corresponding to all parents of i being =1
        Bool=F
        #o is currently a |PA(i)|-length vectors filled with ones, this is one possible assignment to parents of i.
        #we want to cycle through all possible assignments to the parents. Note that each parent can vary
        #from 1 to the domain size of that parent (given in N.Pa)
        #We do this by increasing o lexicographically, i.e. add one to the 'last' entry (closest to the nth
        #entry of o) of o that is still <domain size of the associated parent. For all entries 'after' this
        #one (closer to nth entry), set these values to 1. This is now the next parent assignment lexicographically.
        #Keep doing this until all values of o are equal to the parent domain size (i.e. o=N.Pa), this is the
        #last parent assignment. Moving through the parent assignments in this manner will produce each
        #possible assignment exactly once.
        while(Bool==F){
          W=which(o!=N.Pa)
          if(length(W)==0){
            #If o=N.Pa, we have dealt with all parent assignments and we exit this while loop
            Bool=T
          }
          else{
            #If o=/=N.Pa yet:
            M=max(W)
            if(M==length(PA[[i]])){
              o[length(PA[[i]])]=o[length(PA[[i]])]+1
            }
            else{
              o[M]=o[M]+1
              o[(M+1):length(PA[[i]])]=1
            }
            #Change o to be the next parent assignment lexicographically
          }
          if(Bool==F){
            pref=sample(N[i])
            for(j in 1:N[i]){
              o.f=c(o,j)
              CPT[[i]][matrix(o.f,1)]=pref[j]
            }
            #(If we are not done) for this next parent assignment, o, enter a random ordering of 1-|dom(i)|
            #into CPT(i) as the associated preference order of variable i 
          }
        }
      }
      #CPT(i) is now completely generated and has been added as the next element of the list CPT
    }
    #CPT now contains n complete CPTs (of appropriate sizes for the given CP-net structure, A, and the
    #generated domain sizes, N)
    D=c()
    for(j in 1:n){
      D=c(D,Deg(CPT[[j]],dig))
    }
    #D is a n-length vector. For each CPT (in the list 'CPT') it says whether it is degenerate or not
    if(!any(D=="degenerate")){
      DEGENERATE=FALSE
    }
    #If all of the generated CP-nets are "non-degenerate", we have generated a complete set of valid
    #CPTs and we are done - exit the while loop
    #If any if the CPTs are degenerate, we return to the start of the While loop and generate a new
    #set of CPTs
  }
  return(CPT)
  #We return a full set of (n) valid CPTs for the given CP-net structure (and domain restriction)
}

#rand.cpn
#Input - n, number of variables
#        p, maximum number of parents for any variable
#        d, maximum domain size for any variable
#Output - A random CP-net (with an acyclic sructure) that obeys the input conditions

rand.cpn<-function(n,p,d){
  Pr=c()
  counter=1
  while(length(Pr)<n){
    Next=generate_primes(min=10^(counter-1),max=10^counter)
    Pr=c(Pr,Next)
    counter=counter+1
  }
  Pr=Pr[1:n]
  #Pr - A length n vector of the first n primes in order
  P=prod(Pr^d)
  Bool=F
  i=1
  while(Bool==F){
    if((10^i)>P){
      Bool=T
    }
    i=i+1
  }
  dig=max(c(floor(-log(10^(-i),base=2))+10,100))
  #dig - required level of precision for the functions, given the input values
  L=rand.dag(n,p)
  #L - a random dagcode for a DAG with n nodes, where each node has max p parents
  A=make.dag(L)
  #A - the adjacency matrix of the DAG associated with L
  CPT=rand.CPT(A,d,dig)
  #CPT - a randomly generated list of n valid CPTs for a CP-net with structure A (where no variable
  #has domain larger than d)
  return(list(A=A,CPT=CPT))
  #Output the CP-net with structure A and CPTs given by CPT
  #This is a random CP-net, with an acyclic sructure, that obeys the input conditions
}


#EXAMPLE
#Once the above functions have been loaded, we generate a random CP-Net using rand.cpn
CPN=rand.cpn(4,3,3)
#Generates an acyclic CP-net with four variables, each of which may be either binary or tertiary
#(domain size 2 or 3). Note that my setting p=n-1 we have not set any restriction. Any variable in a
#acyclic 4-variable CP-net can only ever have up to 3 parents.

#An example of the generated CP-Net is given below.
#The first element in the list CPN will be the Adjacency matrix (A) of the CP-net structure
CPN[[1]]
#      [,1] [,2] [,3] [,4]
#[1,]    0    1    1    0
#[2,]    0    0    0    0
#[3,]    0    1    0    0
#[4,]    0    1    1    0

#So, for example, this shows us that variables 1,3, and 4 are all parents of variable 2
#and variables 1 and 4 are parents of variable 3. Variables 1 and 4 have no parents

A=CPN[[1]]

#The second element in CPN is the list of the 4 CPTs. So CPN[[2]][[1]] is the CPT of variable 1,
#CPN[[2]][[2]] is the CPT of variable 2 and so on.
#We will call this list CPT
CPT=CPN[[2]]

CPT[[1]]
#[1] 3 1 2

#As variable 1 has no parents, CPT(1) is simply the preference order over variable 1
#As this list is of length 3, we can see that variable 1 has domain size 3
#Let us say that the domain of variable 1 is {a1,a2,a3}
#The kth entry of the preference order is the preference position of the kth value
#Thus, a1 is the 3rd most preferred, a2 is the most preferred, and a3 is the 2nd most preferred
#So this preference order is a2>a3>a1

#Now, lets say the values of variable 2 are b1,b2,..., variable 3 is c1,c2,.., variable 4 is d1,d2,...
#Recall that we do not yet know which of these variables are binary or tertiary

#We know that variable 2 has variables 1,3,4 as parents. To see the preference order over variable 2
#under the different parental assignments, we use the following command
CPT[[2]][1,1,1,]
#This gives the preference order over variable 2 under
#(Var1=a1,Var3=c1,Var4=d1)
#[1] 2 1

#As the preference order is of length 2, Variable 2 is binary (domain {b1,b2})
#Thus, under this assignment of values to the parents, the preference for variable 2 is b2>b1

dim(CPT[[2]])
#[1] 3 2 2 2

#As we input the parent values to CPT[[2]] in the above way to obtain the preference order,
#looking at the dimension of CPT[[2]] shows us the domain sizes for the parennts (as well as Var2)
#This shows us Var1 is tertiary (as we knew), Var 3&4 are binary, and Var2 is binary (as we knew)

CPT[[2]][3,2,1,]
#[1] 1 2

#So under Var1=a3, Var3=c2, Var4=d1, we have that the preference over Var2 is b1>b2

CPT[[3]][3,2,]
#[1] 2 1

#We know the parents of Var3 are Var1&4 and we know they are tertiary and binary
#This shows us that, under Var1=a3, Var4=d2, the preference over Var3 is c2>c1

#Hopefully, this example of how we access the CPT structures for preference orders, should make 
#the process of the random CPT construction clearer

#In the DominanceTestingFunction.R script we give a function for obtaining the domain sizes of different
#variables (essentially it extracts the last entry in the dimension vector of the CPT array of the var)
#Here, however we know the domain sizes, and we give them in the following vector N

N=c(3,2,2,2)

#An outcome is an assignment of values to all variables. Here, we treat this as an n-vector
#(a length 4 vector for this example). The kth entry will have a value between 1 and the domain size
#of variable k - indicating the value taken by Var k.

o=c(3,1,2,2)

#This is an example of an outcome. In o we have Var1=a3, Var2=b1, Var3=c2, Var4=d2

#The following code randomly generates outcomes for our CP-net, given that N, the vector of domain sizes,
#has been previously calculated

n=4
o1=c()
for(k in 1:n){
  o1[k]=sample(N[k],1)
}

o1
#[1] 3 2 1 2

#This is the way in which we generated random pairs of outcomes for dominance queries.
