classdef TestInterpHermite < TestCase
    %TestInterpHermite Test interpHermite (examples)
    
    properties
        f
        x
        xx
    end
    
    methods
        function self = TestInterpHermite(name)
            self@TestCase(name);
        end
        function setUp(self)
            self.f = @(x) [sin(5*x); 5*cos(5*x); -25*sin(5*x); -125*cos(5*x)];
            self.x = (0:20)/20*2*pi;
            self.xx = (0:360)/360*2*pi;
        end
        function testNoOutput(self)
            Y = self.f(self.x);
            interpHermite(self.x,Y(1,:),Y(2,:),Y(3,:),Y(4,:),'output','sp')
        end
        function testObtainValuesAndDerivatives(self)
            Y = self.f(self.x);
            [vq,vqd,vqdd,vqddd] = interpHermite(self.x,Y(1,:),Y(2,:),Y(3,:),Y(4,:),self.xx);
            assert(norm(vq - sin(5*self.xx)) < 0.2);
        end
        function testVersusPWCH(self)
            Y = self.f(self.x);
            ppInterpHermite = interpHermite(self.x,Y(1,:),Y(2,:),'output','pp');
            ppPWCH = pwch(self.x,Y(1,:),Y(2,:));
            assert(norm(ppInterpHermite.coefs - ppPWCH.coefs) < 1e-10);
        end
        function tearDown(self)
            self.f = [];
            self.x = [];
            self.xx = [];
        end
    end   
end

