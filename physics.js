(function(root){
  'use strict';
  var g=9.81;
  // q=[theta1,theta2,omega1,omega2]; m1=m2=1 kg, L1=L2=1 m.
  // Mass matrix: [[2,cos(d)],[cos(d),1]]. Angles from downward vertical.
  function derivative(q){
    var d=q[0]-q[1],c=Math.cos(d),s=Math.sin(d);
    var r1=-s*q[3]*q[3]-2*g*Math.sin(q[0]);
    var r2=s*q[2]*q[2]-g*Math.sin(q[1]);
    return [q[2],q[3],(r1-c*r2)/(2-c*c),(2*r2-c*r1)/(2-c*c)];
  }
  function add(q,k,h){return q.map(function(v,i){return v+h*k[i];});}
  function step(q,h){
    var a=derivative(q),b=derivative(add(q,a,h/2)),c=derivative(add(q,b,h/2)),d=derivative(add(q,c,h));
    return q.map(function(v,i){return v+h*(a[i]+2*b[i]+2*c[i]+d[i])/6;});
  }
  function energy(q){return q[2]*q[2]+q[3]*q[3]/2+q[2]*q[3]*Math.cos(q[0]-q[1])-g*(2*Math.cos(q[0])+Math.cos(q[1]));}
  function point(q){return [Math.sin(q[0])+Math.sin(q[1]),Math.cos(q[0])+Math.cos(q[1])];}
  function distance(a,b){var x=point(a),y=point(b);return Math.sqrt(Math.pow(x[0]-y[0],2)+Math.pow(x[1]-y[1],2));}
  function integrate(q,h,t){var max=0,e=energy(q);q=q.slice();for(var i=0;i<Math.round(t/h);i++){q=step(q,h);max=Math.max(max,Math.abs(energy(q)-e)/29.43);}return {q:q,drift:max};}
  function verify(){
    var q=[120*Math.PI/180,-10*Math.PI/180,0,0],a=integrate(q,1/240,1),b=integrate(q,1/480,1),c=integrate(q,1/960,1);
    function error(x,y){var n=0;for(var i=0;i<4;i++)n+=Math.pow(x[i]-y[i],2);return Math.sqrt(n);}
    var e1=error(a.q,b.q),e2=error(b.q,c.q),long=integrate(q,1/240,60),equilibrium=step([0,0,0,0],1/240);
    return {drift60:long.drift,error240:e1,error480:e2,ratio:e1/e2,equilibrium:error(equilibrium,[0,0,0,0]),passed:long.drift<0.0001&&e1/e2>10&&e1/e2<22&&error(equilibrium,[0,0,0,0])===0};
  }
  root.Physics={step:step,energy:energy,point:point,distance:distance,verify:verify};
})(this);

