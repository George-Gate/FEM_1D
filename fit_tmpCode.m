kkk=3;
data=[log(nList)/log(10),log(Err{kkk})/log(10)];
data=data(end-5:end,:);
[result,gof]=fit(data(:,1),data(:,2),'poly1');
figure();
plot(result);hold on;
plot(data(:,1),data(:,2),'o');