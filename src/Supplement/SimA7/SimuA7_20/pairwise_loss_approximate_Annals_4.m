function [fun] = pairwise_loss_approximate_Annals_4(n,y,predict_all,w)

% surrogate pairwise loss function 4
fun=0;
for j=1:n
    for i=1:n
        if i~=j    
            m_i=predict_all(i,:);                     %第i个个体的预测值
            y_i=y(i);                                 %第i个个体的真实值
            m_j=predict_all(j,:);                     %第j个个体的预测值
            y_j=y(j);                                 %第j个个体的真实值
            zij=(y_i-y_j)/2;
            if zij==0
                sgn_zij=0;
            elseif zij~=0
                sgn_zij=zij/abs(zij);
            end
            inner_function=m_i-m_j;                        %从内往外第一个函数
            outter_function=-w*(sgn_zij*inner_function)';  %从内往外第二个函数
            if (1+outter_function)>=0
                psi=1+outter_function;
            elseif (1+outter_function)<0
                psi=0;
            end
            fun=fun+min(1,psi);
        elseif i==j
            fun=fun+0;
        end
    end
end

fun=(fun/(n*(n-1)));

end