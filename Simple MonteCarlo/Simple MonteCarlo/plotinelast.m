
hold on
for i = 1:26
   loglog(inelast(:,1),inelast(:,i+1),'DisplayName',num2str(VarName1(i)));
end