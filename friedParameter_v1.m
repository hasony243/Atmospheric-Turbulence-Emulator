A = lambdaWrapped; 
n = N; 
result_distance=zeros(size(A));
vector=A(:);
for r=1:size(A,1)
    for c=1:size(A,2)
        %difference between actual matrix and current position
        %with perfect matches, we need to find 1 or -1 in temp
        temp=abs(A-A(r,c));
        v=cat(3,temp-1,temp+1);
        v=min(abs(v),[],3);
        v(r,c)=inf;%ignore current position as an option
        [~,ind2]=min(v);
        [r2,c2]=ind2sub(size(A),ind2(1));
        result_distance(r,c)=sqrt((r-r2)^2+(c-c2)^2);
    end
end
sample=result_distance(randi(numel(result_distance),1e6,1));
fprintf('the sample mean is %.5f\nthe true mean is   %.5f\n',...
    mean(sample),mean(result_distance(:)))