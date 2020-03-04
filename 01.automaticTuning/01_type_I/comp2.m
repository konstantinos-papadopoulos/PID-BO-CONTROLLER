clc
clear all

% choice = input('pi/pid_3d/pid_2d:','s');
choice = 'pi';
switch lower(choice)
    case 'pi'
        Tp1 = 1;
        Tp2 = 0.01:0.02:1;
        Tp3 = 0.01:0.02:1;
        tp2 = Tp2 / Tp1  ;
        tp3 = Tp3 / Tp1  ;
        Tp2length = length(Tp2);
        Tp3length = length(Tp3);
        kp = 3;
        outerloop = length(Tp2);
        innerloop = length(Tp3);
        w1 = zeros(outerloop,innerloop);
        w2 = zeros(outerloop,innerloop);
        w3 = zeros(outerloop,innerloop);
        numtn = zeros(outerloop,innerloop);
        dentn = zeros(outerloop,innerloop);
        tn = zeros(outerloop,innerloop);
        ti = zeros(outerloop,innerloop);

            for i=1:outerloop
                for j=1:innerloop
                    w1(i,j) = (Tp1 + Tp2(i) + Tp3(j))/Tp1                      ;
                    w2(i,j) = (Tp1*Tp3(j) + Tp1*Tp2(i) + Tp2(i)*Tp3(j))/(Tp1^2);
                    w3(i,j) = Tp1*Tp2(i)*Tp3(j)/(Tp1^3)                        ;
                    numtn(i,j) = w1(i,j)^3 + w3(i,j) - 2*w1(i,j)*w2(i,j)       ;
                    dentn(i,j) = w1(i,j)^2 - w2(i,j)                           ;
                    tn(i,j) = numtn(i,j) /dentn(i,j)                           ;
                    ti(i,j) = 2*kp*(w1(i,j) -tn(i,j) )                         ;
                end
            end
        unity = ones(outerloop,innerloop);
        sizetn = size(tn);
        x_un = ones(innerloop);
        figure(1)
        mesh(tp3,tp2,tn)
        hold on
        mesh(tp3,tp2,unity)
        xlabel('tp3')
        ylabel('tp2')
        zlabel('tn')
    case 'pid_3d'
        Tp1 = 1;
        Tp2 = input('Tp2 = ');
        Tp3 = 0.01:0.01:3;
        Tp4 = 0.01:0.01:2;
        tp1 = Tp1 / Tp1  ;
        tp2 = Tp2 / Tp1  ;
        tp3 = Tp3 / Tp1  ;
        tp4 = Tp4 / Tp1  ;
        kp = 3           ;
        outerloop = length(Tp3);
        innerloop = length(Tp4);
        p1 = zeros(outerloop,innerloop);
        p2 = zeros(outerloop,innerloop);
        p3 = zeros(outerloop,innerloop);
        p4 = zeros(outerloop,innerloop);
        A1 = zeros(outerloop,innerloop);
        B1 = zeros(outerloop,innerloop);
        C1 = zeros(outerloop,innerloop);
        A2 = zeros(outerloop,innerloop);
        B2 = zeros(outerloop,innerloop);
        C2 = zeros(outerloop,innerloop);
        y  = zeros(outerloop,innerloop);
        x  = zeros(outerloop,innerloop);
            for i=1:outerloop
                for j=1:innerloop
                    p1(i,j) = tp1 + tp2 + tp3(i) + tp4(j);
                    p2(i,j) = tp1*tp2 + tp1*tp3(i) + tp1*tp4(j) +...
                              tp2*tp3(i) + tp2*tp4(j) + ...
                              tp3(i)*tp4(j);
                    p3(i,j) = tp1*tp2*tp3(i) + tp1*tp2*tp4(j) + ...
                              tp1*tp3(i)*tp4(j) + ...
                              tp2*tp3(i)*tp4(j) ;
                    p4(i,j) = tp1 * tp2 * tp3(i) * tp4(j);

                    A1(i,j) = p1(i,j)^2 - p2(i,j);
                    B1(i,j) = -p1(i,j)           ;
                    C1(i,j) = p1(i,j)*(p1(i,j)^2 - 2*p2(i,j)) + p3(i,j);
                    A2(i,j) = p2(i,j)^2 - 2*p1(i,j)*p3(i,j) + p4(i,j)  ;
                    B2(i,j) = p3(i,j);
                    C2(i,j) = p1(i,j)*(p2(i,j)^2 - 2*p1(i,j)*p3(i,j) + 2*p4(i,j));

                    y(i,j)  = (A1(i,j)*C2(i,j) - A2(i,j)*C1(i,j))/(A1(i,j)*B2(i,j) - A2(i,j)*B1(i,j));              
                    x(i,j)  = (C1(i,j) - B1(i,j)*y(i,j))/A1(i,j); 
                end
            end
        unity = zeros(outerloop,innerloop)    ;
        unity_ones = ones(outerloop,innerloop);
        x_un = ones(innerloop);
        figure(1)
        mesh(tp4,tp3,y)
        hold on
        mesh(tp4,tp3,unity)
        hold on
        mesh(tp4,tp3,x)
        hold on
        mesh(tp4,tp3,unity_ones)
    otherwise
        

            tp1 = 1;
            tp2 = 0.001:0.001:2;
            tp3 = 0.5 ;
            tp4 = 0.05;

            kp = 1           ;
            outerloop = length(tp2);
            p1 = zeros(outerloop,1);
            p2 = zeros(outerloop,1);
            p3 = zeros(outerloop,1);
            A1 = zeros(outerloop,1)  ;
            B1 = zeros(outerloop,1)  ;
            C1 = zeros(outerloop,1)  ;
            A2 = zeros(outerloop,1)  ;
            B2 = zeros(outerloop,1)  ;
            C2 = zeros(outerloop,1)  ;
            y  = zeros(outerloop,1)  ;
            x  = zeros(outerloop,1)  ;
            D  = zeros(outerloop,1)  ;
            tn = zeros(outerloop,1)  ;
            tv = zeros(outerloop,1)  ;
            ti_pid = zeros(outerloop,1)  ;
            zeroLine = zeros(outerloop,1);
            j = 1;
            
            for i = 1:outerloop
                p1(i) = tp1 + tp2(i) + tp3               ;
                p2(i) = tp1*tp2(i) + tp1*tp3 + tp2(i)*tp3;
                p3(i) = tp1 * tp2(i) * tp3               ;

                A1(i) = p1(i)^2 - p2(i);
                B1(i) = -p1(i)           ;
                C1(i) = p1(i)*(p1(i)^2 - 2*p2(i)) + p3(i);

                A2(i) = p2(i)^2 - 2*p1(i)*p3(i);
                B2(i) = p3(i);
                C2(i) = p1(i)*(p2(i)^2 - 2*p1(i)*p3(i));
                y(i)  = (A1(i)*C2(i) - A2(i)*C1(i))/(A1(i)*B2(i) - A2(i)*B1(i));              
                x(i)  = (C1(i) - B1(i)*y(i))/A1(i);
                D(i) = x(i)^2 - 4*y(i); 
                tn1 = ( x(i) + sqrt(D(i)) ) / 2 ;
                tn2 = ( x(i) - sqrt(D(i)) ) / 2 ;
                 if isreal(tn1) == 0 && isreal(tn2) == 0
                    if tn1 > 0
                       tn(i) = tn1;
                    else
                       tn(i) = tn2;
                    end

                 else
                     tn(i) = tn1;
                 end

                 tv(i) = y(i)/ tn(i) ; 
                
                 if D(i) < 0
                        i_nit(j) = i;
                        j = j + 1   ;
                        y(i)  = y(i);              
                        x(i)  = x(i); 
                 else
                        y(i)  = 0;              
                        x(i)  = 0; 
                 end
%             ti_3(i) = 2*kp*(q1(i) - z1 - x3(i));
            %I Control     
                 
            end
            i = i + 1;
         
            figure(1)
            str = i_nit(1)        ;
            endpl = str + length(i_nit) - 1;
%             plot(tp3(str:endpl),x(str:endpl),'r',tp3(str:endpl),y(str:endpl),'b',tp3,D,'g','LineWidth',2)
            plot(tp2(1:str),tn(1:str),'r--',tp2(1:str),tv(1:str),'k',tp2,D,'g-.','LineWidth',2.5)
%             plot(tp3,ti_3,'r','LineWidth',2)
            current_line_up = zeros(outerloop,1)  ;
            grid
            hold on
            plot(tp2(endpl:outerloop),tn(endpl:outerloop),'r--',tp2(endpl:outerloop),tv(endpl:outerloop),'k',tp2,D,'g-.','LineWidth',2.5)
            hold on          
end