#SERPINZKY 
#EX3
using LinearAlgebra

function serpinzky(start_x,start_y)
    Sx=[];Sy=[];
    append!(Sx,start_x);append!(Sy,start_y);                        
    
    function s(start_x,start_y,Sx,Sy)
        x=[0,0,0,NaN];y=[0,0,0,NaN];
        for i in 1:3
            x[i]=0.5*(start_x[i]+start_x[mod1(i+1,3)])
            y[i]=0.5*(start_y[i]+start_y[mod1(i+1,3)])
        end
        append!(Sx,x);append!(Sy,y);

        if  norm([x[1]-x[2],y[1]-y[2]])<0.005
            return
        end
        t1_x=[x[1],x[2],start_x[2]]; t1_y=[y[1],y[2],start_y[2]];
        t2_x=[x[2],x[3],start_x[3]]; t2_y=[y[2],y[3],start_y[3]];
        t3_x=[x[3],x[1],start_x[1]]; t3_y=[y[3],y[1],start_y[1]];

        s(t1_x,t1_y,Sx,Sy)
        s(t2_x,t2_y,Sx,Sy)
        s(t3_x,t3_y,Sx,Sy)
    end
          
    s(start_x,start_y,Sx,Sy)
    return Sx,Sy
end




#EX 11
#creo funzione piu esterna
function comboss(nn::Int, k::Int)
    S=[];pivot=[];
    N=collect(1:nn)     #vettore [1,.... ,nn]
    recursion(N,k,pivot,S)

    return S
end
#il pivot è il vettore dei numeri che ho fissato. continuo a fissare numeri, e quindi ad allungare il pivot finche
#non arrivo in uno dei casi semplici che so risolvere, ovvero k=1 (guarda add elem) e k=n (una sola combo possibile)
#se non mi trovo in nessuno dei due casi (else) allungo il pivot
function recursion(n::Vector, k::Int, pivot::Vector, S)
    
    #questo è il caso in cui sono arrivato a k=n ho trovato tutte le soluzioni possibili dato il pivot, e accorcio il pivot
    if k==length(n)
        piv=copy(pivot)
        append!(piv,n)
        push!(S,piv)    #aggiungo combo
        pivot=pivot[1:end-1]    #accorcio il pivot
        if length(pivot)==0     #se il pivot è lungo zero ho finito
            return              #FINE
        end
        pivot[end]+=1
        n=collect(pivot[end]+1 : n[end]) #aggiorno la testa del pivot
        recursion(n, k+1, pivot,S)      #rilancio la funzione con pivot piu corto e quindi con k+1
   
    elseif k==1
        add_elements(S,pivot,n)
        pivot[end]+=1   #update pivot
        recursion(n[2:end], k, pivot,S) 
                        
    else
        pivot_temp=copy(pivot)
        append!(pivot_temp,n[1])
        recursion(n[2:end], k-1, pivot_temp,S)
    end
end

#questa funzione serve a creare le combo e ad aggiungerle al vettore soluzione S 
#una volta che siamo arrivati alla situazione (n k=1) ove tutti gli n-1 elementi sono fissati v1 e faccio combo con
#i restanti v2
function add_elements(S,v1::Vector,v2::Vector)
    for i in v2
        v1_temp=copy(v1)
        append!(v1_temp,i)
        push!(S,v1_temp)
    end
end

###################################################################

#EX 15
#the idea is to find the biggest square, count 1 square operation and than to use the 
#function on both the number that is squared and the rest
function power_efficiency(n::Int)
    res=Dict()
    res["square"]=0
    res["molt"]=0
    res["2"]=false #it's a flag that once it is true, num(2)=2 when 0 num(1)=1   try f(4) to unders
    function num(n::Int, res::Dict)
        if n==1
            res["molt"]+=1 #caso noto
            return
        elseif n==0 
            return
        elseif n==2 #il 2 è tricky, costa 1 operazione la prima volta e 2 per le altre
            if res["2"]==false
                res["2"]=true
                res["square"]+=1
                return 
            else
                res["molt"]+=2 #it could be also 1molt and 1 square, its the same
                return 
            end
        else  
            i=1
            while i*i<=n
                i+=1
            end
            i=i-1 #because it does an extra cycle
            r= n-i*i
            if i != 1
                res["square"]+=1 #being 1**2 =1 there is no squaring
            end
            num(i, res)
            num(r, res)
            return res
        end
    end
    return num(n, res)
end
##################################################################################
#EX 16
#minimizzare la differenza tra sx e dx 
struct knaps
    diff::Float16                   
    v::Vector
    w::Vector
end

function knapsack(v::Vector)

    w=[]; dif=Inf;
    function knap(v,w,dif)

        cost=[]
        for i in 1:length(v)
            v1=copy(v); w1=copy(w)
            append!(w1,v1[i]);  deleteat!(v1,i)

            dif_1=abs(sum(v1)-sum(w1))

            if dif_1==0
                return knaps(dif_1,v1,w1)
            elseif dif_1>dif
                return knaps(dif,v,w)
            end

            push!(cost, knap(v1,w1,dif_1) )
        end
        
        index=argmin( [cost[i].diff for i in 1:length(cost)] )
        return cost[index] 

    end
    ciao=knap(v,w,dif)
    return ciao
end

sum!