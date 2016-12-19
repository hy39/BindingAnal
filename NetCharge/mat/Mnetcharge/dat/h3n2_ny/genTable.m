Ch = Viruses
find(isnan(Ch(:,1)))
Ch(find(isnan(Ch(:,1))),:) = [];
Ch(find(Ch(:,1)==0),:) = [];
%agegroupmin = [0, 10, 20, 40, 65, 80];
%agegroupmax = [10,20, 40, 65, 80, 110];
agegroupmin = [0, 10, 20, 30, 40, 50, 60, 70, 80];
agegroupmax = [10,20, 30, 45, 50, 60, 70, 80,110];
%agegroupmin = [0, 20, 40, 60, 80];
%agegroupmax = [20,40, 60, 80,100];
for i=1:length(agegroupmin)
   if i == 0 
      disp('imposible')
   end
   Ch(find(Ch(:,1)>agegroupmin(i) & Ch(:,1)<=agegroupmax(i)),10)=i;
end

csvwrite('Chlist9.dat', Ch);

% calculate the percentage of the high netcharge
med_charge = 17;
ch_age = [];
sizeage_each = [];
for x=1:length(agegroupmin)
    sizeage_each(x) = length(find(Ch(:,10)>=x & Ch(:,10)<(x+1)));
    ch_age(x,1) = length(find(Ch(:,10)>=x & Ch(:,10)<(x+1) & Ch(:,3) > med_charge)) /  length(find(Ch(:,10)>=x & Ch(:,10)<(x+1)));

    pbin = ch_age(x,1);
    sizeage(x) = sizeage_each(x)
    y0 = pbin;
    options = optimoptions('fsolve','Display','off'); % Option to display output
    ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
    ch_age(x,2) = lb;
    ch_age(x,3) = ub;
end

sizeage
figure
hold on
bar(1:length(agegroupmin), ch_age(:,1));
%errorbar(1:length(agegroupmin),mean_velocity,std_velocity,'.')


% ---- below is not used

%pbin = 0.6226;
%sizeage = sizeage_each(5);
%y0 = 0.4;
%options = optimoptions('fsolve','Display','off'); % Option to display output
%ub = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.025','x','pbin','sizeage'),y0,[options],pbin,sizeage);
%lb = fsolve(inline('binocdf(round(pbin*sizeage),sizeage,x)-0.975','x','pbin','sizeage'),y0,[options],pbin,sizeage);
