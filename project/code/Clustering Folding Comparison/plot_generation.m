%variables vs. epsilon, americasProblem
figure(1) 
hold on
variables=[24 20 20 16 14 14;24 23 22 22 23 22;24 24 24 24 24 24];
D=[2;3;4];
epsilon=[0 0.3 0.4 0.5 0.6 0.7];
s1=scatter(epsilon,variables(1,:)./variables(1,1));
s2=scatter(epsilon,variables(2,:)./variables(2,1));
s3=scatter(epsilon,variables(3,:)./variables(3,1));
legend([s1 s2 s3],'2-coloring','3-coloring','4-coloring')
xlim([0,1])
ylim([0,1])
xlabel('epsilon in epsilon net')
ylabel('variables (% of original problem)')