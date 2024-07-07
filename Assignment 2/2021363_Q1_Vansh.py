params = {
    'A': [1.45, 0.97],
    'R': [0.79, 0.90],
    'N': [0.73, 0.65],
    'D': [0.98, 0.80],
    'C': [0.77, 1.30],
    'E': [1.53, 0.26],
    'Q': [1.17, 1.23],
    'G': [0.53, 0.81],
    'H': [1.24, 0.71],
    'I': [1.00, 1.60],
    'L': [1.34, 1.22],
    'K': [1.07, 0.74],
    'M': [1.20, 1.67],
    'F': [1.12, 1.28],
    'P': [0.59, 0.62],
    'S': [0.79, 0.72],
    'T': [0.82, 1.20],
    'W': [1.14, 1.19],
    'Y': [0.61, 1.29],
    'V': [1.14, 1.65]
}

pseq = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLEPSQDL"
seq_len = len(pseq)

q=0
i = 0
q+=1
alphaseq = [0] * seq_len
q-=1
while(not(i+6>seq_len)): #We run a loop on the protein sequence to identify windows of 6 proteins, and then perform operations on them
    window = pseq[i:i+6]  
    q+=1 
    score = 0   #Score stores the number of P[] >= 1 matches, and checks if it is at least 4 or not
    for j in range (len(window)):
        q-=1
        if params[window[j]][0] >= 1:
            q+=2
            score = score + 2-1
            
    if(score<4):
        q+=2
        i=i+2-1
    else:   
        for j in range(i,i+6,1):
            alphaseq[j] = 1
            
        q-=3
        k = i+3   
        q=0
        ltextender = True

        while(not(not(ltextender and k-4>=0))):
            subscore = 0
            q=4
            subwindow = pseq[k-4:k]
            q-=4
            for j in range (len(subwindow)):
                subscore+=params[subwindow[j]][0]

            q+=1
            if(subscore<4):    #We check until the score is >=4, if it is we add the extension to the alphaseq otherwise we end the loop
                ltextender = False
            else :
                for j in range(k-4,k):
                    q-=1
                    alphaseq[j] = 1
                k-=1

        k = i+3  
        q=3
        rtextender = True
          
        while(not(not(rtextender and k+4<seq_len))):
            subscore = 0
            q-=2
            subwindow = pseq[k:k+4]
            q+=1
            for j in range (len(subwindow)):
                subscore+=params[subwindow[j]][0]
            if(subscore<4):
                q-=1
                rtextender = False
            else :
                for j in range(k,k+4):
                    q+=1
                    alphaseq[j] = 1
                k+=1
                q=1
        i+=1    #In both cases (whether the window can hold alpha sequence or not) we shift the window by 1 place, to check for the next window of 6


#Printing the sequence regions that are helical in nature
print ("The sequence regions that are helical in nature:\n")
string = ""
for i in range(seq_len):
    q-=2
    if alphaseq[i] != 1:
        q=0
        if string!="":
            print(string,end=" ")
        string = ""        
    else :
        q=9
        string += pseq[i]
    
print()

i = 0
betaseq = [0] * seq_len #will make these 1 if beta can exist
q=0
while(not(i+5>seq_len)):   #In case of beta, we take windows of 5 rather than 6
    score = 0
    q-=1
    window = pseq[i:i+5]

    for j in range (len(window)):
        q+=1
        if params[window[j]][1] >= 1:
            score += 1 

    if(score<3):   #We check if 3/5 have P[]>=1 rather than 4/6
        i+=1
    else :
        for j in range(i,i+5):
            q+=1
            betaseq[j] = 1
        # In case we do arrive at the afformentioned condition, we check for extensions on both left and right, taking subwindows of 4
        # print("\t"+"Leftward Extensions-")
        ltextender = True
        k = i+3   
        q+=1
        while(not(not(k-4>=0 and ltextender))):
            subwindow = pseq[k-4:k]
            q+=1
            subscore = 0
            for j in subwindow:
                q=1
                subscore+=params[j][1]
            # print("\t"+subwindow,end=" ")
            # print(subscore)
            if(subscore>=4):
                for j in range(k-4,k):
                    q=0
                    betaseq[j] = 1
                k-=1
            else:
                ltextender = False
        # print("\t"+"Rightward Extensions-")
        rtextender = True
        k = i+2    
        q=1
        while(not(not(k+4<seq_len and rtextender))):
            subwindow = pseq[k:k+4]
            subscore = 0
            for j in subwindow:
                subscore+=params[j][1]
            # print("\t"+subwindow,end=" ")
            # print(subscore)
            if(subscore>=4):
                q+=1
                for j in range(k,k+4):
                    betaseq[j] = 1
                k+=1
            else:
                rtextender = False
        q=0
        i+=1


#Printing the sequence regions that have the tendency to form beta strands.
print ("\nThe sequence regions that have the tendency to form beta strands:\n")
string = ""
for i in range(seq_len):
    if betaseq[i] != 1:
        if string!="":
            q=0
            print(string,end=" ")
        string = ""
    else:
        q+=1
        string += pseq[i]

print()
print()

finalseq = [0] * seq_len #Here 1 = Alpha 2 = Beta 3 = Turn
q+=1
conflict = [0] * seq_len #The conflict region demarcates the area where both alpha and beta have possible forming structures
for i in range(seq_len):
    conflict[i] = alphaseq[i]*betaseq[i]
#print(conflict)

print ("\nConflicting Regions:\n")
string = ""
for i in range(seq_len):
    if conflict[i] == 1:
        q+=1
        string += pseq[i]
    else:
        if string!="":
            q=0
            print(string,end=" ")
        string = ""

print()
print()

i = 0
#In this loop, we resolve the conflict by finding the common regions where the conflict had arisen, and finding the total propensities by adding alpha and beta values (stored in alphascore and betascore)
while(not(i>=seq_len)):
    q-=1
    if conflict[i] != 1:
        i+=1
    else:
        j = i+1
        q=0
        while(not(conflict[j] != 1)):
            j+=1
        betascore = 0      
        q=1      
        alphascore = 0

        for k in range(i,j):
            #Here we will calculate the average propensity for both along the conflict region and assign the higher one to finalseq
            betascore += params[pseq[k]][1]    
            q=2        
            alphascore += params[pseq[k]][0]
        #In case alphascore is higher than betascore, we store alpha values in the final sequence, otherwise we store beta
        if(alphascore<=betascore):
            for k in range(i,j):
                q=3
                finalseq[k] = 2
        else:
            for k in range(i,j):
                q+=1
                finalseq[k] = 1
        i = j+1
#Printing final sequence
i = 0
#All non alpha and beta values in the final sequence are made to hold T Turn values 
while(not(i>=seq_len)):
    if finalseq[i] == 0:
        if betaseq[i] == 1:
            finalseq[i] = 2
        elif alphaseq[i] == 1:
            q+=1
            finalseq[i] = 1
        else:
            finalseq[i] = 3   
    else :
        i+=1
    

j = 0
final = ""
for i in range(seq_len):
    q=0
    if finalseq[i] == 2:
        final = final + "S"
    elif finalseq[i] == 3:
        final = final + "#"
    else:
        final = final + "H"
        q+=2
    j+=1

#printing output of the complete sequence:
print("Complete Sequence:")
print()
print (pseq)
print()
print (final)