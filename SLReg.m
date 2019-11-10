classdef SLReg < handle
    %SLReg Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LAMBDA
        M
        Mp
        H
        f
        J
        psi
        N
        Gz
        D
        Nu
        u_prev
        yzad
        deltaup
        umin
        umax
        ymin
        ymax
    end
    
    methods
        function obj = SLReg(D, N, Nu, lambda, Umin, Umax, Ymin, Ymax)
            %DM Construct an instance of this class
            %   Detailed explanation goes here
            obj.N = N;
            obj.D = D;
            obj.Nu = Nu;
            obj.psi = eye(N);
            obj.LAMBDA = eye(Nu)*lambda;
            obj.u_prev = 73;
            obj.umin = ones(Nu,1).*Umin;
            obj.umax = ones(Nu,1).*Umax;
            obj.J = tril(Nu);
            obj.ymin = ones(N,1).*Ymin;
            obj.ymax = ones(N,1).*Ymax;
        end
        
        function [] = reset(obj,u_p)
            obj.u_prev = u_p;
            obj.deltaup=zeros(1,obj.D-1)';
        end
        
        function [] = setValue(obj, yzad)
            yzad(1:obj.N)= yzad;
            obj.yzad=yzad';
        end
        
        function [] = newModel(obj, FDp)
            h = getLinearModel(obj.u_prev, FDp);
            A = zeros(2,2);
            B = zeros(2,2);
            C = zeros(1,2);
            D1 = zeros(1,2);
            
            C1 = 0.35;  
            C2 = 0.3;
            alfa1 = 20;
            alfa2 = 22;
            tau = 150;
            
            A(1,1) = alfa1 * h(1)^(-5/2)/(2*C1) - 2/(3*C1) * (obj.u_prev + FDp) * h(1)^(-3);
            A(2,1) = alfa1/(6*C2) * h(1)^(-1/2) * h(2)^-2;
            A(2,2) = -2 * alfa1/(3*C2) * sqrt(h(1)) * h(2)^-3 + alfa2/(2*C2) * h(2)^(-5/2);
            B(1,1) = h(1)^(-2)/(3*C1);
            B(1,2) = h(1)^(-2)/(3*C1);
            C(1,2) = 1;
            sys = ss(A,B,C,D1,'InputDelay',[tau,0]);
            
            obj.Gz = c2d(tf(sys), 0.5, 'zoh');
            
            
            s = step(obj.Gz, obj.D/obj.Gz.Ts);
            s = s(:,:,1);
            
            % Wyznaczanie macierzy M
            obj.M = zeros(obj.N, obj.Nu);
            for i=1:1:obj.Nu
                obj.M(i:obj.N,i)=s(1:obj.N-i+1)';
            end


            % Wyznaczanie macierzy Mp
            obj.Mp= zeros(obj.N, obj.D-1);
            for i=1:obj.N
               for j=1:obj.D-1
                  if i+j<=obj.D
                     obj.Mp(i,j)=s(i+j)-s(j);
                  else
                     obj.Mp(i,j)=s(obj.D)-s(j);
                  end   
               end
            end

            % Wektor wzmocnieñ
            obj.H=2*(obj.M'*obj.psi*obj.M+obj.LAMBDA);
        end
        
        function u = countValue(obj,y)
            %count actual value of output.
            %   Detailed explanation goes here
            
            % aktualizacja wektora aktualnej wartoœci wyjœcia
            yk=ones(obj.N,1)*y;

            % wyliczenie nowego wektora odpowiedzi swobodnej
            y0=yk+obj.Mp*obj.deltaup;
    
            obj.f = -2*obj.M'*obj.psi*(obj.yzad-y0);
            
            % wyliczenie wektora zmian sterowania
            deltauk(1) = fmincon(@(x) 1/2*x'*obj.H*x+obj.f'*x, obj.u_prev, [ -obj.J; obj.J; -obj.M ; obj.M], [-obj.umin+obj.deltaup(1) ; obj.umax-obj.deltaup(1) ; vertcat(-obj.ymin +y0, obj.ymax- y0)], [], [], obj.umin, obj.umax);

            % prawo regulacji
            u= obj.u_prev + deltauk(1);
            obj.u_prev = u;

            % aktualizacja poprzednich zmian sterowania
            obj.deltaup=[deltauk(1) obj.deltaup(1:end-1)']';
        end
    end
end

