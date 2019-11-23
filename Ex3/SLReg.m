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
        MFs
        s
        D
        Nu
        size
        u_prev
        yzad
        deltaup
        umin
        umax
        ymin
        ymax
        MembershipFunctions
    end
    
    methods
        function obj = SLReg(D, N, Nu, lambda, Umin, Umax, Ymin, Ymax, Gzs, MembershipFunctions)
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
            obj.MembershipFunctions = MembershipFunctions;
            obj.size = size(Gzs);
            obj.s = cell(obj.size);
            %Creates steps responses for all transfer functions.
            for i=1:obj.size(1,1)
                for j=1:obj.size(1,2)
                    s = step(Gzs{i,j}, D/Gzs{i,j}.Ts);
                    s = s(:,:,1);
                    obj.s{i,j} = s(1:D);
                end
            end
        end
        
        function [] = reset(obj,u_p)
            obj.u_prev = u_p;
            obj.deltaup=zeros(1,obj.D-1)';
        end
        
        function [] = setValue(obj, yzad)
            yzad(1:obj.N)= yzad;
            obj.yzad=yzad';
        end
        
        function u = countValue(obj,y)
            %count actual value of output.
            %   Detailed explanation goes here

            num = zeros(1,obj.D);
            den = zeros(1,obj.D);

            % Count step response.
            for i=1:obj.size(1,1)
               for j= 1:obj.size(1,2)
                    v = obj.MembershipFunctions{i,j}.getValue(y);
                    num = num + v.*obj.s{i,j};
                    den = den + v;
               end
            end
            s = num./den;

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

            % Wektor wzmocnie�
            obj.H=2*(obj.M'*obj.psi*obj.M+obj.LAMBDA);

            % aktualizacja wektora aktualnej warto�ci wyj�cia
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

