colors={'r','g','b'};
refcolors={'--r','--g','--b'};
v=zeros(2,size(vectors,2));

figure(1)
hold on
for i=1:size(vectors,1)
    v(2,:)=vectors(i,:);
    %plot(v(:,1),v(:,2),colors{j(i)})
    plot3(v(:,1),v(:,2),v(:,3),colors{j(i)})
end
for i=1:size(x,1)
    v(2,:)=x(i,:)/norm(x(i,:));
    %plot(v(:,1),v(:,2),refcolors{i})
    plot3(v(:,1),v(:,2),v(:,3),refcolors{i})
end
title('Clustering: color denotes variable assignment, dashed lines are centers')
view(3)
saveas(1,'coloring_n6_clustering.png')

figure(2)
hold on
for i=1:size(vectors,1)
    v(2,:)=vectors(i,:);
    %plot(v(:,1),v(:,2),colors{j(i)})
    plot3(v(:,1),v(:,2),v(:,3),colors{J(i)})
end
for i=1:size(centers,1)
    v(2,:)=centers(i,:);
    %plot(v(:,1),v(:,2),refcolors{i})
    plot3(v(:,1),v(:,2),v(:,3),refcolors{i})
end
title('Classifying: color denotes variable assignment, dashed lines are centers')
view(3)
saveas(2,'coloring_n6_classifying.png')