function [x,y,MI_test] = MI_heatmap_search(MI)


for iInits = 1:10000

       x = randi(size(MI,1));
       y = randi(size(MI,2));

       n_runs = 1000;

       MI_test = MI(x,y);


       for iRun = n_runs
           x_test = 0;
           y_test = 0;
           try
               if MI(x-1,y) > MI_test
                   MI_test = MI(x-1,y);
                   x_test = x-1;
                   y_test = y;
               end
           end

           try
               if MI(x-1,y-1) > MI_test
                   MI_test = MI(x-1,y-1);
                   x_test = x-1;
                   y_test = y-1;
               end
           end

           try
               if MI(x-1,y+1) > MI_test
                   MI_test = MI(x-1,y+1);
                   x_test = x-1;
                   y_test = y+1;
               end
           end

           try
               if MI(x,y+1) > MI_test
                   MI_test = MI(x,y+1);
                   x_test = x;
                   y_test = y+1;
               end
           end

           try
               if MI(x,y-1) > MI_test
                   MI_test = MI(x,y-1);
                   x_test = x;
                   y_test = y-1;
               end
           end

           try
               if MI(x+1,y+1) > MI_test
                   MI_test = MI(x+1,y+1);
                   x_test = x+1;
                   y_test = y+1;
               end
           end

           try
               if MI(x+1,y) > MI_test
                   MI_test = MI(x+1,y);
                   x_test = x+1;
                   y_test = y;
               end
           end

           try
               if MI(x+1,y-1) > MI_test
                   MI_test = MI(x+1,y-1);
                   x_test = x+1;
                   y_test = y-1;
               end
           end       

           if ~x_test
                x = x_test;
                y = y_test;
           end

       end
       
       a(iInits,:) = [x y];
end

a(find(a(:,1)==0 & a(:,2)==0),:) = [];

[u,~,c] = unique(a,'rows');
[~,ix]=max(accumarray(c,1));
mdrows=u(ix,:);

end

