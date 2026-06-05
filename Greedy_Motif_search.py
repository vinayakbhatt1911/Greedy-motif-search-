#Greedy Motif Search Algorithm 
#Position weight matrix Function
def build_profile(motif):
   k=len(motif[0])
   pwm=[]
   for i in range(k):
      count={'A':1,'T':1,'G':1,'C':1}
      for m in motif:
         count[m[i]]+=1
      total=sum(count.values())
      prob={}
      for x in 'ATGC':
       prob[x]=count[x]/total
      pwm.append(prob)
   return  pwm
#Score calculate Function
def  calculate_score(motif):
   score=0
   for i in range(len(motif[0])):
      A = T = G = C = 0
    
      for m in motif:
        if(m[i]=='A'):
           A+=1
        if(m[i]=='T'):
           T+=1
        if(m[i]=='G'):
           G+=1
        if(m[i]=="C"):
           C+=1
      max_count=max(A,T,G,C)
      score=score+(len(motif)-max_count)

   return score                   
# Find best K-mer Function
def profile_most_probable_kmer(seq,k,profile):
   best=seq[0:k]
   max_p=-1
   for i in range(len(seq)-k+1):
      kmer=seq[i:i+k]
      p=calculate_probability(kmer,profile)

      if(p>max_p):
         max_p=p
         best=kmer
    
   return best
#Get probability Function
def calculate_probability(kmer,profile):
   
   p=1
   for i in range(len(kmer)):
      ch=kmer[i]
      p=p*profile[i][ch]
    
   return p

# (GREEDY FUNCTION)
def greedy(DNA,k):
   best=[]

   for seq in DNA:
      best.append(seq[0:k])

   first=DNA[0]

   for i in range(len(first)-k+1):
      motif=[]
      motif.append(first[i:i+k])

      for j in range(1,len(DNA)):
         profile=build_profile(motif)
         m=profile_most_probable_kmer(DNA[j],k,profile)
         motif.append(m)

      if calculate_score(motif) < calculate_score(best):
         best=list(motif)
    
   return best
      

# MAIN FUNCTION 
n=int(input("Enter how many numbers of DNA sequence are there:"))
DNA=[]
for i in range(n):
   seq=input("Enter DNA sequence:")
   DNA.append(seq)
k=int(input("Enter k value:"))

result=greedy(DNA,k)
print("The Best Motif is:")
print("Motif:",result)
print("Score",calculate_score(result))
