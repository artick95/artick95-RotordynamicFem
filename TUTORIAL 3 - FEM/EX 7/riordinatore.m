function[A]=riordinatore(A,O,ind)
%in O, su una stessa riga ci stanno gli indici delle righe e delle colonne
%da scambiare: ad ex prima riga di O: [a b]========> scambia riga a con
%riga b e colonna a con colonna b

%ind: se ind=1 scambia sia righe che colonna, altrimenti scambia solo le
%righe

for k=1:size(O,1)
    a=O(k,1);
    b=O(k,2);
      
    r=A(a,:);
    A(a,:)=A(b,:);
    A(b,:)=r;
    
    if ind==1
        c=A(:,a);
        A(:,a)=A(:,b);
        A(:,b)=c;
    end
end