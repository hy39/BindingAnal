figure;
Mnew = M;

for i = 1:33
  if M(i,1) < 0
        M(i,5) = -1;
  end
  if M(i,1) == 0
        M(i,5) = 0;
  end
  if M(i,1) > 0
        M(i,5) = 1;
  end
end

for i = 1:33
  for j = i:33
   
        if i~=j
		if Mnew(i,1) == Mnew(j,1)
	           if abs(Mnew(i,4)-Mnew(j,4))<0.2
                      Mnew(j,1) = Mnew(j,1)+0.05
                   end
		end	 
        end
  end
end
csvwrite('Mnewlist.dat',Mnew);
csvwrite('Mlist.dat',M);
%scatter(Mnew(:,1),Mnew(:,4),'filled','SizeData',16);
%hold on;

idx1 = find(M(:,1)==-1);
x1 = Mnew(idx1,4);

idx2 = find(M(:,1)==0);
x2 = Mnew(idx2,4);

idx3 = find(M(:,1)==1);
x3 = Mnew(idx3,4);

grp = [zeros(1,length(x1)),ones(1,length(x2)),2*ones(1,length(x3))];

boxplot([x1',x2',x3'],grp);
hold on;

plot(Mnew(:,1),Mnew(:,4),'.');