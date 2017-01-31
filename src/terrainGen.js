//Dual contouring code based on implementation by Nick Gildea (which was based on a number of other sources)
//All other code written by Lane Spangler (las4vc)


//Constants
var QEF_SWEEPS = 4;
var QEF_ERROR = .000001;

var MATERIAL_AIR = 0;
var MATERIAL_DIRT = 1;
var MATERIAL_ROCK = 2;
var MATERIAL_GRASS = 3;

var CHUNK_SIZE = 64;

var CUBOID_CACHE_LIMIT = 64;
var SPHERE_CACHE_LIMIT = 64;
var NOISE_CACHE_LIMIT = 256;

var MAX_INTERSECTIONS = 6;
var HOMOGENOUSLEAF = 0;
var LEAF = 1;
var INTERIOR = 2;

//Inheritance
function extend(base, sub) {

  var origProto = sub.prototype;
  sub.prototype = Object.create(base.prototype);
  for (var key in origProto)  {
     sub.prototype[key] = origProto[key];
  }
  Object.defineProperty(sub.prototype, 'constructor', { 
    enumerable: false, 
    value: sub 
  });
}


function printV(v){
	console.log(v.x+", "+v.y+", "+v.z);
}


//Matrices etc for solving qef
function fnorm(m){
	var a = m;

	if(!a.m10){
		a = new Mat3();
		a.copySymmetric(m);
	}

	return Math.sqrt((a.m00*a.m00) + (a.m01*a.m01) + (a.m02*a.m02)
		+ (a.m10*a.m10) + (a.m11*a.m11) + (a.m12 * a.m12)
		+ (a.m20*a.m20) + (a.m21*a.m21) * (a.m22*a.m22));
}

function off(m){
	var a = m;
	if(!a.m10){
		return Math.sqrt(2*((a.m01*a.m01) + (a.m02 * a.m02) + (a.m12 * a.m12)));
	}

	return Math.sqrt((a.m01*a.m01) + (a.m02 * a.m02) + (a.m10 * a.m10)
		+ (a.m12*a.m12) + (a.m20*a.m20) + (a.m21 * a.m21));

}

function mmul(a,b){
	var m = new Mat3();
	m.set(a.m00 * b.m00 + a.m01 * b.m10 + a.m02 * b.m20,
		a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21,
		a.m00 * b.m02 + a.m01 * b.m12 + a.m02 * b.m22,
		a.m10 * b.m00 + a.m11 * b.m10 + a.m12 * b.m20,
		a.m10 * b.m01 + a.m11 * b.m11 + a.m12 * b.m21,
		a.m10 * b.m02 + a.m11 * b.m12 + a.m12 * b.m22,
		a.m20 * b.m00 + a.m21 * b.m10 + a.m22 * b.m20,
		a.m20 * b.m01 + a.m21 * b.m11 + a.m22 * b.m21,
		a.m20 * b.m02 + a.m21 * b.m12 + a.m22 * b.m22);
	return m;
}

function mmul_ata(a){
	var m = new SMat3();
	m.setSymmetric(a.m00 * a.m00 + a.m10 * a.m10 + a.m20 * a.m20,
	     a.m00 * a.m01 + a.m10 * a.m11 + a.m20 * a.m21,
	     a.m00 * a.m02 + a.m10 * a.m12 + a.m20 * a.m22,
	     a.m01 * a.m01 + a.m11 * a.m11 + a.m21 * a.m21,
	     a.m01 * a.m02 + a.m11 * a.m12 + a.m21 * a.m22,
	     a.m02 * a.m02 + a.m12 * a.m12 + a.m22 * a.m22);
	return m;
}

function transpose(a){
	var m = new Mat3();
	m.set(a.m00, a.m10, a.m20, a.m01, a.m11, a.m21, a.m02, a.m12, a.m22);
	return m;
}

function vmul(a, v){
	return new THREE.Vector3((a.m00 * v.x) +(a.m01*v.y) + (a.m02*v.z), 
			(a.m10*v.x) + (a.m11*v.y) + (a.m12*v.z),
			(a.m20*v.x) + (a.m21*v.y) + (a.m22*v.z));
}

function vmul_sym(m,v){
	var a = new Mat3();
	a.copySymmetric(m);

	return vmul(a,v);
}

function givensRot01(m, c, s){
    m.set(c * m.m00 - s * m.m01, s * m.m00 + c * m.m01, m.m02, c * m.m10 - s * m.m11,
          s * m.m10 + c * m.m11, m.m12, c * m.m20 - s * m.m21, s * m.m20 + c * m.m21, m.m22);
}

function givensRot02(m, c, s){
    m.set(c * m.m00 - s * m.m02, m.m01, s * m.m00 + c * m.m02, c * m.m10 - s * m.m12, m.m11,
          s * m.m10 + c * m.m12, c * m.m20 - s * m.m22, m.m21, s * m.m20 + c * m.m22);
}

function givensRot12(m, c, s){
    m.set(m.m00, c * m.m01 - s * m.m02, s * m.m01 + c * m.m02, m.m10, c * m.m11 - s * m.m12,
          s * m.m11 + c * m.m12, m.m20, c * m.m21 - s * m.m22, s * m.m21 + c * m.m22);
}

function givensCoefficients(a_pp, a_pq, a_qq,co){
		co.c = 0;
		co.s = 0;

		if(a_pq == 0){
			co.c = 1;
			return co;
		}

		var tau = (a_qq - a_pp) / (2 * a_pq);
		var stt = Math.sqrt(1.0 + tau*tau);
		var tan = 1.0 / ((tau >= 0) ? (tau + stt) : (tau - stt));
		co.c = 1.0 / Math.sqrt(1 + tan*tan);
		co.s = tan * co.c;
		return co;
}

function schurRot01(m,co){
	givensCoefficients(m.m00,m.m01,m.m11,co);

	var cc = co.c*co.c;
	var ss = co.s*co.s;
	var mix = 2 * co.c * co.s * m.m01;
	m.setSymmetric(cc * m.m00 - mix + ss * m.m11, 0, co.c * m.m02 - co.s * m.m12,
                       ss * m.m00 + mix + cc * m.m11, co.s * m.m02 + co.c * m.m12, m.m22);
}

function schurRot02(m,co){
	givensCoefficients(m.m00,m.m02,m.m22,co);

	var cc = co.c*co.c;
	var ss = co.s*co.s;
	var mix = 2 * co.c * co.s * m.m02;
    m.setSymmetric(cc * m.m00 - mix + ss * m.m22, co.c * m.m01 - co.s * m.m12, 0,
                   m.m11, co.s * m.m01 + co.c * m.m12, ss * m.m00 + mix + cc * m.m22);
}

function schurRot12(m,co){
	givensCoefficients(m.m11,m.m12,m.m22,co);

	var cc = co.c*co.c;
	var ss = co.s*co.s;
	var mix = 2 * co.c * co.s * m.m12;
    m.setSymmetric(m.m00, co.c * m.m01 - co.s * m.m02, co.s * m.m01 + co.c * m.m02,
                   cc * m.m11 - mix + ss * m.m22, 0, ss * m.m11 + mix + cc * m.m22);
}

function rotate01(vtav, v){
	if(vtav.m01 == 0){
		return;
	}

	var co = {};
	schurRot01(vtav,co);
	givensRot01(v,co.c,co.s);
}

function rotate02(vtav, v){
	if(vtav.m02 == 0){
		return;
	}

	var co = {};
	schurRot02(vtav,co);
	givensRot02(v,co.c,co.s);
}

function rotate12(vtav, v){
	if(vtav.m12 == 0){
		return;
	}

	var co = {};
	schurRot12(vtav,co);
	givensRot12(v,co.c,co.s);
}

function getSymmetricSvd(a, vtav, v, tol, max_sweeps){
	vtav.copySymmetric(a);
	v.set(1,0,0,
		0,1,0,
		0,0,1);
	var delta = tol*fnorm(vtav);
	var i;
	for(i = 0; i < max_sweeps && off(vtav) > delta; i++){
		rotate01(vtav,v);
		rotate02(vtav,v);
		rotate12(vtav,v);
	}
}

function calcError(A, x, b){
	var vtmp = vmul(A,x);
	vtmp = vtmp.subVectors(b,vtmp);
	return vtmp.dot(vtmp);
}

function symCalcError(origA, x, b){
	var A = new Mat3();
	A.copySymmetric(origA);
	var vtmp = vmul(A,x);
	vtmp = vtmp.subVectors(b,vtmp);
	return vtmp.dot(vtmp);
}

function pinv(x, tol){
	return (Math.abs(x) < tol || Math.abs(1.0 / x) < tol) ? 0 : (1.0 / x);
}

function pseudoinverse(d, v, tol){
	var m = new Mat3();
	var d0 = pinv(d.m00, tol);
	var d1 = pinv(d.m11, tol);
	var d2 = pinv(d.m22, tol);
	m.set(v.m00 * d0 * v.m00 + v.m01 * d1 * v.m01 + v.m02 * d2 * v.m02,
	    v.m00 * d0 * v.m10 + v.m01 * d1 * v.m11 + v.m02 * d2 * v.m12,
	    v.m00 * d0 * v.m20 + v.m01 * d1 * v.m21 + v.m02 * d2 * v.m22,
	    v.m10 * d0 * v.m00 + v.m11 * d1 * v.m01 + v.m12 * d2 * v.m02,
	    v.m10 * d0 * v.m10 + v.m11 * d1 * v.m11 + v.m12 * d2 * v.m12,
	    v.m10 * d0 * v.m20 + v.m11 * d1 * v.m21 + v.m12 * d2 * v.m22,
	    v.m20 * d0 * v.m00 + v.m21 * d1 * v.m01 + v.m22 * d2 * v.m02,
	    v.m20 * d0 * v.m10 + v.m21 * d1 * v.m11 + v.m22 * d2 * v.m12,
	    v.m20 * d0 * v.m20 + v.m21 * d1 * v.m21 + v.m22 * d2 * v.m22);
	return m;
}

function solveSymmetric(A, b, x, svd_tol, svd_sweeps, pinv_tol){
	var mtmp = new Mat3();
	var pinv = new Mat3();
	var V = new Mat3();
	var VTAV = new SMat3();
	getSymmetricSvd(A,VTAV, V, svd_tol, svd_sweeps);
	pinv = pseudoinverse(VTAV, V, pinv_tol);
	var r = vmul(pinv,b);
	x.x = r.x;
	x.y = r.y;
	x.z = r.z;
	return symCalcError(A, x, b);
}

function solveLeaseSquares(a, b, x, svd_tol, svd_sweeps, pinv_tol){
	var at = new Mat3();
	var ata = new SMat3();
	var atb = new THREE.Vector3(0,0,0);
	at = transpose(a);
	ata = mmul_ata(a);
	atb = vmul(at,b);
	return solveSymmetric(ata, atb, x, svd_tol, svd_sweeps, pinv_tol);
}

SMat3 = function(){
	this.setSymmetric(0,0,0,0,0,0);
};

SMat3.prototype.copySymmetric = function(sym){
	this.setSymmetric(sym.m00,sym.m01,sym.m02,sym.m11,sym.m12,sym.m22);
};

SMat3.prototype.setSymmetric = function(m00,m01,m02,m11,m12,m22){
	this.m00 = m00;
	this.m01 = m01;
	this.m02 = m02;
	this.m11 = m11;
	this.m12 = m12;
	this.m22 = m22;
};

Mat3 = function(){
	this.set(0,0,0,0,0,0,0,0,0);
};

Mat3.prototype.set = function(m00,m01,m02,m10,m11,m12,m20,m21,m22){
	this.m00 = m00;
	this.m01 = m01;
	this.m02 = m02;
	this.m10 = m10;
	this.m11 = m11;
	this.m12 = m12;
	this.m20 = m20;
	this.m21 = m21;
	this.m22 = m22;
};

Mat3.prototype.copySymmetric = function(sym){
	this.setSymmetric(sym.m00,sym.m01,sym.m02,sym.m11,sym.m12,sym.m22);
};

Mat3.prototype.setSymmetric = function(a00,a01,a02,a11,a12,a22){
	this.set(a00, a01, a02, a01, a11, a12, a02, a12, a22);
};

//Quadratic error function for calculating vertex positions and determining if nodes can be collapsed

QEF = function(){
	this.reset();
};

QEF.prototype.getMassPoint = function(){
	return this.massPoint;
};

QEF.prototype.add  = function(qef){
	this.ata_00 += qef.ata_00;
	this.ata_01 += qef.ata_01;
	this.ata_02 += qef.ata_02;
	this.ata_11 += qef.ata_11;
	this.ata_12 += qef.ata_12;
	this.ata_22 += qef.ata_22;
	this.atb_x += qef.atb_x;
	this.atb_y += qef.atb_y;
	this.atb_z += qef.atb_z;
	this.btb += qef.btb;
	this.massPoint_x += qef.massPoint_x;
	this.massPoint_y += qef.massPoint_y;
	this.massPoint_z += qef.massPoint_z;
	this.numPoints += qef.numPoints;
};

QEF.prototype.setAta = function(){
	this.ata = new SMat3();
	this.ata.setSymmetric(this.ata_00,this.ata_01,this.ata_02,
		this.ata_11,this.ata_12,this.ata_22);
};

QEF.prototype.setAtb = function(){
	this.atb = new THREE.Vector3(this.atb_x,this.atb_y,this.atb_z);
};

QEF.prototype.solve = function(outx, svd_tol, svd_sweeps, pinv_tol){
	if(this.numPoints == 0){
		console.log("Can't solve with no intersections");
	}

	this.massPoint = new THREE.Vector3(this.massPoint_x,this.massPoint_y,this.massPoint_z);
	this.massPoint.multiplyScalar(1.0/this.numPoints);

	this.setAta();
	this.setAtb();

	var tmpv = vmul_sym(this.ata,this.massPoint);
	this.atb.sub(tmpv);

	this.x = new THREE.Vector3(0,0,0);

	var result = solveSymmetric(this.ata,this.atb,this.x,svd_tol,svd_sweeps,pinv_tol);
	this.x.add(this.massPoint);
	this.setAtb();

	this.hasSolution = true;

	outx.x = this.x.x;
	outx.y = this.x.y;
	outx.z = this.x.z;

	return result;

};

QEF.prototype.initialize = function(intersectionCount, normals, positions){
	this.hasSolution = false;
	var i;
	for(i = 0; i < 12; i++){

		if(normals[i]){
			var n = normals[i];
			var p = positions[i];
			this.ata_00 += n.x*n.x;
			this.ata_01 += n.x*n.y;
			this.ata_02 += n.x*n.z;
			this.ata_11 += n.y*n.y;
			this.ata_12 += n.y*n.z;
			this.ata_22 += n.z*n.z;
			var dot = n.dot(p);
			this.atb_x += dot*n.x;
			this.atb_y += dot*n.y;
			this.atb_z += dot*n.z;
			this.btb += dot*dot;
			this.massPoint_x += p.x;
			this.massPoint_y += p.y;
			this.massPoint_z += p.z;
			this.numPoints++;
		}
	}
};

QEF.prototype.getError = function(){
	if(!this.hasSolution){
		console.log("Can't find error if there's no solution");
	}

	return this.getPosError(this.x);
};

QEF.prototype.getPosError = function(p){
	if(!this.hasSolution){
		this.setAta();
		this.setAtb();
	}

	var atax = vmul_sym(this.ata,p);
	return p.dot(atax) - 2*p.dot(this.atb) + this.btb;
};

QEF.prototype.reset = function(){
	this.hasSolution = false;
	this.ata_00 = 0;
	this.ata_01 = 0;
	this.ata_02 = 0;
	this.ata_11 = 0;
	this.ata_12 = 0;
	this.ata_22 = 0;
	this.atb_x = 0;
	this.atb_y = 0;
	this.atb_z = 0;
	this.btb = 0;
	this.massPoint_x = 0;
	this.massPoint_y = 0;
	this.massPoint_z = 0;
	this.numPoints = 0;
};


//Data needed for heterogenous nodes (the only nodes that are actually drawn)
NodeDrawInfo = function(){
	this.index = -1;
	this.corners = [0,0,0,0,0,0,0,0];
	this.position = new THREE.Vector3(0,0,0);
	this.normal = new THREE.Vector3(0,0,0);
	this.positions = [];
	this.normals = [];
	this.qef = new QEF();
};

//Octree data structure
OctreeNode = function(){
	this.type = -1; //Types are None, Interior, HomogenousLeaf, Leaf
	this.min = new THREE.Vector3(0,0,0); //The position of the cube corner with the smallest x,y,z values
	this.size = 0; //Node dimensions are size x size x size
	this.drawInfo = null;
	this.children = []; //8 children in an octree
};

OctreeNode.prototype.expandLeaf = function(){
}

//Position offsets for each child node
var childOffsets = [new THREE.Vector3(0,0,0),
					new THREE.Vector3(0,0,1),
					new THREE.Vector3(0,1,0),
					new THREE.Vector3(0,1,1),
					new THREE.Vector3(1,0,0),
					new THREE.Vector3(1,0,1),
					new THREE.Vector3(1,1,0),
					new THREE.Vector3(1,1,1)];

var edgevmap = [
	[0,4],[1,5],[2,6],[3,7],
	[0,2],[1,3],[4,6],[5,7],
	[0,1],[2,3],[4,5],[6,7]
];

Vertex = function(pos,normal){
	this.xyz = pos;
	this.normal = normal;
};


Chunk = function(min, hermiteData){
	this.root = new OctreeNode();
	this.root.size = CHUNK_SIZE;
	this.root.min = min;
	this.root.type = INTERIOR;

	this.df = hermiteData;
}

Chunk.prototype.drawOctree = function(node, scene){

	if(!node){
		return;
	}

	if(node.type == LEAF){
		var material = new THREE.MeshBasicMaterial({color: 0x00ff00});
		material.wireframe = true;

		var geometry = new THREE.BoxGeometry(node.size,node.size,node.size);

		var mesh = new THREE.Mesh(geometry,material);
		scene.add(mesh);
		mesh.position.set(node.min.x+node.size/2,node.min.y+node.size/2,node.min.z+node.size/2);

	}

	var i;
	for(i = 0; i < 8; i++){
		this.drawOctree(node.children[i],scene);
	}
};

//To build a leaf node, we need to sample the density function at each of the node corners
//If there are any sign changes, we need to compute intersections, normals, and minimize the QEF to get the vertex
Chunk.prototype.buildLeaf = function(node){
	var corners = [];
	var materialChange = false;
	var emptyCorner = false;
	var currentMaterial = 0;
	var i;

	for(i = 0; i < 8; i++){
		var cornerPos = new THREE.Vector3();
		cornerPos.addVectors(node.min,new THREE.Vector3(childOffsets[i].x*node.size,childOffsets[i].y*node.size,childOffsets[i].z*node.size));
		var material = this.df.getMaterialData(cornerPos.x,cornerPos.y,cornerPos.z);

		if(material == MATERIAL_AIR){
			emptyCorner = true;
		}

		if(i == 0){
			currentMaterial = material;
		}else if(material != currentMaterial){
			materialChange = true;
		}
		corners[i] = material;
	}

	if(!materialChange){
		//Homogenous node, does need intersection data
		node.type = HOMOGENOUSLEAF;
		node.material = material;
		return node;
	}

	var max_intersections = MAX_INTERSECTIONS;
	var positions = [];
	var normals = [];
	var intersectionCount = 0;

	for(i = 0; i < 12 && intersectionCount < max_intersections; i++){
		
		var c1 = edgevmap[i][0];
		var c2 = edgevmap[i][1];


		var m1 = corners[c1];
		var m2 = corners[c2];

		if (m1 == m2)
		{
			// no zero crossing on this edge
			continue;
		}

		//Otherwise find the intersection point
		var p1 = new THREE.Vector3();
		p1.addVectors(node.min,childOffsets[c1]);
		var p2 = new THREE.Vector3();
		p2.addVectors(node.min,childOffsets[c2]);
		var intersection = this.df.getSurfaceIntersection(p1,p2);

		//Store position of intersection and normal at intersection for each edge
		positions[i] = intersection;
		normals[i] = this.df.getSurfaceNormal(intersection);

		intersectionCount++;
	}
	var drawInfo = new NodeDrawInfo();
	drawInfo.qef.initialize(intersectionCount, normals, positions);
	drawInfo.position = new THREE.Vector3(0,0,0);
	drawInfo.qef.solve(drawInfo.position,QEF_ERROR, QEF_SWEEPS, QEF_ERROR);

	var min = node.min;
	var max = new THREE.Vector3();
	max.addVectors(min, new THREE.Vector3(node.size,node.size,node.size));

	if(drawInfo.position.x < min.x || drawInfo.position.x > max.x ||
		drawInfo.position.y < min.y || drawInfo.position.y > max.y ||
		drawInfo.position.z < min.z || drawInfo.position.z > max.z){
		drawInfo.position = drawInfo.qef.getMassPoint();
	}

	for(i = 0; i < 12; i++){
		if(normals[i]){
			drawInfo.normal.add(normals[i]);
		}
	}

	drawInfo.normal.normalize();

	drawInfo.corners = corners;
	drawInfo.positions = positions;
	drawInfo.normals = normals;

	node.type = LEAF;
	node.drawInfo = drawInfo;

	return node;
};

//Top down octree construction
//Note: it would be probably be faster (but less simple) to do this construction from the leaf nodes and up
//Generate leaf nodes, then generate parent nodes breadth first until node size is CHUNK_SIZE
Chunk.prototype.buildOctree = function(node, ignoreBounds){
	if(!node){
		return null;
	}

	if(node.size == this.df.leafSize){
		return this.buildLeaf(node);
	}

	var childSize = node.size / 2; //The children of a 4x4x4 cube are 2x2x2 etc.
	//console.log("Child size " + childSize);
	var heterogenous = false;
	var i;

	
	for(i = 0; i < 8; i++){
		var child = new OctreeNode();
		child.size = childSize;
		child.min = new THREE.Vector3();
		var offset = new THREE.Vector3(childOffsets[i].x*childSize,childOffsets[i].y*childSize,childOffsets[i].z*childSize);
		child.min.addVectors(node.min,offset);
		child.type = INTERIOR;

		//If the operation is unbounded then build out the octree depth first
		if(!this.df.AABB || ignoreBounds){
			node.children[i] = this.buildOctree(child,true); //Recursively build out the children of the child
		}else{
			//Otherwise we only continue on nodes that intersect the operation's bounds
			var childAABB = new AABB(child.min,new THREE.Vector3(childSize,childSize,childSize));

			if(this.df.AABB.contains(childAABB)){
				//If the child is completely contained in the operation bounds then there's no need to keep checking its children
				node.children[i] = this.buildOctree(child,true);
			}
			if(childAABB.intersects(this.df.AABB)){
				//Otherwise the recursive call happens only when the child intersects the operation bounds
				node.children[i] = this.buildOctree(child,false);
			}else{
				child.type = HOMOGENOUSLEAF;
				child.material = MATERIAL_AIR;
				node.children[i] = child;
			}
		}



		heterogenous = heterogenous | (node.children[i].type != HOMOGENOUSLEAF);
	}
		


	if(!heterogenous){
		node.type = HOMOGENOUSLEAF;
		node.material = node.children[0].material;
		node.children = [];
	}

	return node;
};

Chunk.prototype.simplifyOctree = function(node, threshold){

	if(node.type != INTERIOR){
		return node;
	}

	var qef = new QEF();
	var signs = [-1,-1,-1,-1,-1,-1,-1,-1];
	var midsign = -1;
	var canCollapse = true;
	//var edgeCount = 0;

	var i;
	for(i = 0; i < 8; i++){
		node.children[i] = this.simplifyOctree(node.children[i],threshold);
		var child = node.children[i];
		if(child.type != HOMOGENOUSLEAF){

			if(child.type == INTERIOR){
				canCollapse = false;
			}else{
					qef.add(child.drawInfo.qef);
					midsign = child.drawInfo.corners[7-i];
					signs[i] = child.drawInfo.corners[i];

					//edgeCount++;
			}
		}else{
			signs[i] = child.material;
			midsign = child.material;
		}
	}

	//Wait to return until all of the children have been simplified
	if(!canCollapse){
		return node;
	}

	var qefPosition = new THREE.Vector3(0,0,0);
	qef.solve(qefPosition,QEF_ERROR,QEF_SWEEPS,QEF_ERROR);
	var error = qef.getError();

	if(error > threshold){
		return node;
	}

	if(qefPosition.x < node.min.x || qefPosition.x > (node.min.x + node.size) ||
		qefPosition.y < node.min.y || qefPosition.y > (node.min.y + node.size) ||
		qefPosition.z < node.min.z || qefPosition.z > (node.min.z + node.size)){
		var mp = qef.getMassPoint();
		qefPosition = new THREE.Vector3(mp.x,mp.y,mp.z);
	}

	var drawInfo = new NodeDrawInfo();

	for(i = 0; i < 8; i++){
		if(signs[i] == -1){
			drawInfo.corners[i] = midsign;
		}else{
			drawInfo.corners[i] = signs[i];
		}
	}

	drawInfo.normal = new THREE.Vector3(0,0,0);
	for(i =0; i < 8; i++){
		if(node.children[i]){
			var child = node.children[i];
			if(child.type == LEAF){
				drawInfo.normal.add(child.drawInfo.normal);

				//extrapolate intersection normals and positions
				var j;
				for(j = 0; j < 12; j++){	
					var childEdge = child.drawInfo.normals[j];
					var edge = drawInfo.normals[j];

					if(childEdge){
						if(!edge){
							drawInfo.normals[j] = child.drawInfo.normals[j];
							drawInfo.positions[j] = child.drawInfo.positions[j];
						}else{
							drawInfo.normals[j].add(child.drawInfo.normals[j]);
							drawInfo.normals[j].normalize();
							drawInfo.positions[j].add(child.drawInfo.positions[j]);
						}
					}
				}
			}
		}
	}
	drawInfo.normal.normalize();
	drawInfo.position = qefPosition;
	drawInfo.qef = qef;


	for(i = 0; i < 8; i++){
		node.children[i] = null;
	}

	node.type = LEAF;
	node.drawInfo = drawInfo;

	return node;
};

function copyOctree(node, copy){
	copy.min = new THREE.Vector3(node.min.x,node.min.y,node.min.z);
	copy.size = node.size;
	copy.type = node.type;

	if(node.type == INTERIOR){
		var i;
		for(i = 0; i < 8; i++){
			copy.children[i] = new OctreeNode();
			copy.children[i] = this.copyOctree(node.children[i],copy.children[i]);
		}
	}else if(node.type == HOMOGENOUSLEAF){
		copy.material = node.material;
		copy.children = null;
	}else if(node.type == LEAF){
		var drawInfo = new NodeDrawInfo();

		drawInfo.qef.add(node.drawInfo.qef);

		drawInfo.position = new THREE.Vector3(node.drawInfo.position.x,node.drawInfo.position.y,node.drawInfo.position.z);
		drawInfo.normal = new THREE.Vector3(node.drawInfo.normal.x,node.drawInfo.normal.y,node.drawInfo.normal.z);

		drawInfo.index = node.drawInfo.index;

		var i;
		for(i = 0; i < 8; i++){
			drawInfo.corners[i] = node.drawInfo.corners[i];
		}

		for(i = 0; i < 12; i++){
			if(node.drawInfo.normals[i]){
				drawInfo.normals[i] = new THREE.Vector3(node.drawInfo.normals[i].x,node.drawInfo.normals[i].y,node.drawInfo.normals[i].z);
				drawInfo.positions[i] = new THREE.Vector3(node.drawInfo.positions[i].x,node.drawInfo.positions[i].y,node.drawInfo.positions[i].z);
			}
		}

		copy.drawInfo = drawInfo;
	}

	return copy;
}

Chunk.prototype.generateSimplifiedMesh = function(simplificationThreshold){
	if(!this.root){
		return;
	}

	this.vertexBuffer = [];
	this.indexBuffer = [];
	this.materialBuffer = [];
	this.vertexCount = 0;
	this.indexCount = 0;
	this.faceCount = 0;

	//copy root into simpleRoot
	this.simpleRoot = new OctreeNode();
	copyOctree(this.root,this.simpleRoot);
	console.log(JSON.stringify(this.simpleRoot));

	this.simpleRoot = this.simplifyOctree(this.simpleRoot,simplificationThreshold);

	this.generateVertexIndices(this.simpleRoot);
	this.contourCellProc(this.simpleRoot);
	this.simpleRoot = null;
};

Chunk.prototype.generateVertexIndices = function(node){
	if(!node){
		return;
	}

	if(node.type == INTERIOR){
		var i;
		for(i = 0; i < 8; i++){
			this.generateVertexIndices(node.children[i]);
		}
	}

	if(node.type == LEAF){
		if(!node.drawInfo){
			console.log("No draw info for node");
		}
		node.drawInfo.index = this.vertexCount;
		var v = new Vertex(node.drawInfo.position,node.drawInfo.normal);

		var i;
		this.vertexBuffer[this.vertexCount] = v;
		this.vertexCount++;

	}
};

var cellProcFaceMask = [[0,4,0],[1,5,0],[2,6,0],[3,7,0],[0,2,1],[4,6,1],[1,3,1],[5,7,1],[0,1,2],[2,3,2],[4,5,2],[6,7,2]];
var cellProcEdgeMask = [[0,1,2,3,0],[4,5,6,7,0],[0,4,1,5,1],[2,6,3,7,1],[0,2,4,6,2],[1,3,5,7,2]];

Chunk.prototype.contourCellProc = function(node){
	if(node.type == HOMOGENOUSLEAF){
		return;
	}

	if(node.type == INTERIOR){
		var i,j;
		for(i = 0; i < 8; i++){
			this.contourCellProc(node.children[i]);
		}

		for(i = 0; i < 12; i++){
			var faceNodes = [];
			var c = [cellProcFaceMask[i][0], cellProcFaceMask[i][1]];

			faceNodes[0] = node.children[c[0]];
			faceNodes[1] = node.children[c[1]];

			this.contourFaceProc(faceNodes,cellProcFaceMask[i][2]);
		}

		for(i = 0; i < 6; i++){
			var edgeNodes = [];
			var c = [cellProcEdgeMask[i][0],
				cellProcEdgeMask[i][1],
				cellProcEdgeMask[i][2],
				cellProcEdgeMask[i][3]];

			for(j = 0; j < 4; j++){
				edgeNodes[j] = node.children[c[j]];
			}

			this.contourEdgeProc(edgeNodes, cellProcEdgeMask[i][4]);
		}
	}
};

var faceProcFaceMask = [[[4,0,0],[5,1,0],[6,2,0],[7,3,0]],
						[[2,0,1],[6,4,1],[3,1,1],[7,5,1]],
						[[1,0,2],[3,2,2],[5,4,2],[7,6,2]]];

var faceProcEdgeMask = [[[1,4,0,5,1,1],[1,6,2,7,3,1],[0,4,6,0,2,2],[0,5,7,1,3,2]],
						[[0,2,3,0,1,0],[0,6,7,4,5,0],[1,2,0,6,4,2],[1,3,1,7,5,2]],
						[[1,1,0,3,2,0],[1,5,4,7,6,0],[0,1,5,0,4,1],[0,3,7,2,6,1]]];

var orders = [[0,0,1,1],[0,1,0,1]];
Chunk.prototype.contourFaceProc = function(node,dir){
	if(node[0].type == HOMOGENOUSLEAF || node[1].type == HOMOGENOUSLEAF){
		return;
	}

	if(node[0].type == INTERIOR ||
		node[1].type == INTERIOR){

		var i,j;
		for(i = 0; i < 4; i++){
			var faceNodes = [];
			var c = [faceProcFaceMask[dir][i][0],faceProcFaceMask[dir][i][1]];

			for(j = 0; j < 2; j++){
				if(node[j].type != INTERIOR){
					faceNodes[j] = node[j];
				}else{
					faceNodes[j] = node[j].children[c[j]];
				}
			}

			this.contourFaceProc(faceNodes, faceProcFaceMask[dir][i][2]);
		}

		for(i = 0; i < 4; i++){
			var edgeNodes = [];
			var c = [faceProcEdgeMask[dir][i][1],
				faceProcEdgeMask[dir][i][2],
				faceProcEdgeMask[dir][i][3],
				faceProcEdgeMask[dir][i][4]];

			var order = orders[faceProcEdgeMask[dir][i][0]];

			for(j = 0; j < 4; j++){
				if(node[order[j]].type == LEAF){
					edgeNodes[j] = node[order[j]];
				}else{
					edgeNodes[j] = node[order[j]].children[c[j]];
				}
			}

			this.contourEdgeProc(edgeNodes, faceProcEdgeMask[dir][i][5]);
		}

	}
};

var edgeProcEdgeMask = [[[3,2,1,0,0],[7,6,5,4,0]],
						[[5,1,4,0,1],[7,3,6,2,1]],
						[[6,4,2,0,2],[7,5,3,1,2]]];

Chunk.prototype.contourEdgeProc =function(node, dir){
	if(node[0].type == HOMOGENOUSLEAF || node[1].type == HOMOGENOUSLEAF || node[2].type == HOMOGENOUSLEAF || node[3].type == HOMOGENOUSLEAF){
		return;
	}

	if(node[0].type != INTERIOR &&
		node[1].type != INTERIOR &&
		node[2].type != INTERIOR &&
		node[3].type != INTERIOR){
		this.contourProcessEdge(node,dir);
	}else{
		var i,j;
		for(i = 0; i < 2; i++){
			var edgeNodes = [];
			var c = [edgeProcEdgeMask[dir][i][0],
			edgeProcEdgeMask[dir][i][1],
			edgeProcEdgeMask[dir][i][2],
			edgeProcEdgeMask[dir][i][3]];

			for(j = 0; j < 4; j++){
				if(node[j].type == LEAF){
					edgeNodes[j] = node[j];
				}else{
					edgeNodes[j] = node[j].children[c[j]];
				}
			}

			this.contourEdgeProc(edgeNodes, edgeProcEdgeMask[dir][i][4]);
		}

	}
};

var processEdgeMask = [[3,2,1,0],[7,5,6,4],[11,10,9,8]];
//Actually connect 4 vertices to make a quad
Chunk.prototype.contourProcessEdge = function(node, dir){
	var minSize = 1000000;
	var minIndex = 0;
	var indices = [-1,-1,-1,-1];
	var flip = false;
	var signChange = [false,false,false,false];
	var material = [1,1,1,1];

	var minM1 = 0;
	var minM2 = 0;

	var i;
	for(i = 0; i < 4; i++){
		var edge = processEdgeMask[dir][i];
		var c1 = edgevmap[edge][0];
		var c2 = edgevmap[edge][1];
		var m1 = node[i].drawInfo.corners[c1];
		var m2 = node[i].drawInfo.corners[c2];
		//TODO (probably move to shader) material should be a blend of the two materials
		if(m1 != MATERIAL_AIR){
			material[i] = m1;
		}else{
			material[i] = m2;
		}


		if(node[i].size < minSize){
			minSize = node[i].size;
			minIndex = i;
			minM1 = m1;
			minM2 = m2;
			//if m1 is solid then the quad needs to face the opposite direction
			flip = m1 != MATERIAL_AIR;
		}

		indices[i] = node[i].drawInfo.index;

		signChange[i] = (m1 == MATERIAL_AIR && m2 != MATERIAL_AIR) || (m1 != MATERIAL_AIR && m2 == MATERIAL_AIR);

	}

	if(signChange[minIndex]){

		//console.log("M1: "+minM1+", M2: "+minM2);
		//For now find dominating material or just pick a non air material
		var tri1mat,tri2mat;
		if((material[0] == material[1] || material[0] == material[3]) && material[0] != MATERIAL_AIR ){
			tri1mat = material[0];
		}else if(material[1] == material[3] && material[1] != MATERIAL_AIR){
			tri1mat = material[1];
		}else{
			var i;
			for(i = 0; i < 4; i++){
				if(material[i] != MATERIAL_AIR){
					tri1mat = material[i];
					break;
				}
				tri1mat = MATERIAL_ROCK;
			}
		}

		if((material[0] == material[2] || material[0] == material[3]) && material[0] != MATERIAL_AIR){
			tri2mat = material[0];
		}else if(material[2] == material[3] && material[2] != MATERIAL_AIR){
			tri2mat = material[2];
		}else{
			var i;
			for(i = 0; i < 4; i++){
				if(material[i] != MATERIAL_AIR){
					tri2mat = material[i];
					break;
				}
				tri2mat = MATERIAL_ROCK;
			}
		}

		this.materialBuffer[this.faceCount] = tri1mat;
		this.faceCount++;
		this.materialBuffer[this.faceCount] = tri2mat;
		this.faceCount++;

		if(tri1mat == MATERIAL_AIR || tri2mat == MATERIAL_AIR){
			console.log("Undefined material index");
		}

		//TODO (shader) face should be a blend of each vertex material
		if(!flip){
			this.indexBuffer[this.indexCount] = indices[0];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[1];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[3];
			this.indexCount++;

			this.indexBuffer[this.indexCount] = indices[0];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[3];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[2];
			this.indexCount++;
		}else{
			this.indexBuffer[this.indexCount] = indices[0];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[3];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[1];
			this.indexCount++;

			this.indexBuffer[this.indexCount] = indices[0];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[2];
			this.indexCount++;
			this.indexBuffer[this.indexCount] = indices[3];
			this.indexCount++;
		}
	}
};

//create CSGNodes, give them a hierarchy and operations,
//ask for the result of the root node
//It should return world data

//Let a world be a list of axis aligned, uniform, chunks with coordinates

//Thus low level operations create worlds

//Mid level operations combine worlds

//High level operations combine operations

CSGNode = function(params){
	this.children = [];
	this.childCount = 0;
	this.params = params;
	this.world = null;
}

//Abstract
CSGNode.prototype.operation = function(a,b){
	return this;
}

CSGNode.prototype.addChild = function(child){
	this.children[this.childCount] = child;
	this.childCount++;
}

//Get the output world of the current CSG tree
CSGNode.prototype.getResult = function(){

	if(this.childCount == 0){
		//If a node was not given children then the operation should generate a world
		return this.operation();
		//And the result is just the node with the world field filled out
	}
	var i;


	if(this.childCount == 1){
			return this.operation(this.children[0],null);
	}

	var result = this.operation(this.children[0],this.children[1]);

	for(i = 2; i < this.childCount; i++){
		result = this.operation(result, this.children[i]);
	}

	return result;
}

//Concat -> Add a as an operand to b
CSGNode.prototype.concat = function(a,b){
	a.getResult();

	if(!b){
		return a;
	}

	//get the result of a and add it as a child of b
	b.addChild(a);

	//Return the result of b
	b.getResult();
	return	b;
}

//Combination operations

//Union -> Combine volumes, left operand overwrites right operand
UnionNode = function(params){
	CSGNode.call(this,params);
}

UnionNode.prototype.operation = function(a,b){
	var u = new CSGNode();

	var unionWorld;
	unionWorld = new Chunk(new THREE.Vector3(0,0,0),null);

	if(!a && b){
		unionWorld.root = b.getResult().world.root;
		u.world = unionWorld;
		return u;
	}else if(a && !b){
		unionWorld.root = a.getResult().world.root;
		u.world = unionWorld;
		return u;
	}else if(!a && !b){
		return null;
	}


	var dominantWorld = a.getResult().world;
	var otherWorld = b.getResult().world;

	unionWorld.root = this.getUnion(dominantWorld.root,otherWorld.root);

	u.world = unionWorld;
	return u;
}

//Recursive function to actually find the union of two octrees filled with hermite data
//Assumes trees start at the same depth... could be modified to just let the larger sized nodes dominate until the sizes are equal
UnionNode.prototype.getUnion = function(a,b){

	if(a.type == INTERIOR && b.type == INTERIOR){
		var i;
		for(i = 0; i < 8; i++){
			a.children[i] = this.getUnion(a.children[i],b.children[i]);
		}
			
		return a;
	}

	if(a.type == INTERIOR && b.type == LEAF){
		b.expandLeaf();
		//return this.getUnion(a,b)
		return a;
	}

	if(b.type == INTERIOR && a.type == LEAF){
		a.expandLeaf();
		//return this.getUnion(a,b)
		return a;
	}

	if(a.type == HOMOGENOUSLEAF){
		if(a.material == MATERIAL_AIR){
			return b;
		}else{
			//If the dominant node is non-air and homogenous then it replace any heterogenous node
			return a;
		}
	}

	if(a.type == LEAF && b.type == HOMOGENOUSLEAF){
		//If b is a homogenous node and a is a leaf then we replace any empty corners in A with the material of b
		if(b.material != MATERIAL_AIR){
			var i;
			for(i = 0; i < 8; i++){
				if(a.drawInfo.corners[i] == MATERIAL_AIR){
					a.drawInfo.corners[i] = b.material;
				}
			}
		}
		return a;
	}


	if(a.type == INTERIOR && b.type == HOMOGENOUSLEAF){
		//If b is air then just return a
		if(b.material == MATERIAL_AIR){
			return a;
		}else{
			//Otherwise we recursively get the union of a with the homogenous node
			var i;
			for(i = 0; i < 8; i++){
				a.children[i] = this.getUnion(a.children[i],b);
			}
			return a;
		}
	}

	//Both nodes are leaves then replace the empty corners of a with material b and empty edges of A with edges of B
	if(a.type == LEAF && b.type == LEAF){
		var i;

		for(i = 0; i < 8; i++){


			if(a.drawInfo.corners[i] == MATERIAL_AIR){
				a.drawInfo.corners[i] = b.drawInfo.corners[i];
			}
		}

		//Should there be an intersection limit on a union?
		var qefNeedsUpdate = false;
		var intersectionCount = 0;
		for(i = 0; i < 12; i++){
			var edgeA = a.drawInfo.normals[i];
			var edgeB = b.drawInfo.normals[i];
			if(!edgeA && edgeB){
				a.drawInfo.positions[i] = b.drawInfo.positions[i];
				a.drawInfo.normals[i] = b.drawInfo.normals[i];
				a.drawInfo.normal.add(b.drawInfo.normals[i]);
				intersectionCount++;
				qefNeedsUpdate = true;
			}else if(edgeA){
				intersectionCount++;
			}
			
		}

		if(qefNeedsUpdate){
			a.drawInfo.qef.initialize(intersectionCount,a.drawInfo.normals,a.drawInfo.positions);
			var qefPosition = new THREE.Vector3();
			a.drawInfo.qef.solve(qefPosition,QEF_ERROR, QEF_SWEEPS, QEF_ERROR);

			if(qefPosition.x < a.min.x || qefPosition.x > (a.min.x + a.size) ||
			qefPosition.y < a.min.y || qefPosition.y > (a.min.y + a.size) ||
			qefPosition.z < a.min.z || qefPosition.z > (a.min.z + a.size)){
				var mp = a.drawInfo.qef.getMassPoint();
				qefPosition = new THREE.Vector3(mp.x,mp.y,mp.z);
			}

			a.drawInfo.position = qefPosition;
			a.drawInfo.normal.normalize();
		}

		return a;
	}

};

extend(CSGNode,UnionNode);

DiffNode = function(params){
	CSGNode.call(this,params);
}
//Diff -> Take a volume chunk out of the right operand where the left operand is
DiffNode.prototype.operation = function(a,b){
	var d = new CSGNode();

	var diffWorld;
	diffWorld = new Chunk(new THREE.Vector3(0,0,0),null);

	if(!a && b){
		diffWorld.root = b.getResult().world.root;
		d.world = diffWorld;
		return d;
	}else if(a && !b){
		//Should be air
		diffWorld.root = a.getResult().world.root;
		d.world = diffWorld;
		return d;
	}else if(!a && !b){
		return null;
	}


	var worldToSubtract = a.getResult().world;
	var dominantWorld = b.getResult().world;

	diffWorld.root = this.getDifference(worldToSubtract.root,dominantWorld.root);

	d.world = diffWorld;
	return d;
}

//b-a
DiffNode.prototype.getDifference = function(a,b){
	if(a.type == INTERIOR && b.type == INTERIOR){
		var i;
		for(i = 0; i < 8; i++){
			b.children[i] = this.getDifference(a.children[i],b.children[i]);
		}

		return b;
	}

	if(a.type == INTERIOR && b.type == LEAF){
		console.log("Diff interior and leaf");
		b.expandLeaf();
		return b;
	}

	if(b.type == INTERIOR && a.type == LEAF){
		console.log("Diff interior and leaf");
		a.expandLeaf();
		return b;
	}

	if(a.type == HOMOGENOUSLEAF){
		if(a.material != MATERIAL_AIR){
			a.material = MATERIAL_AIR;
			return a;
		}else{
			return b;
		}
	}

	if(b.type == HOMOGENOUSLEAF){
		if(b.material != MATERIAL_AIR){
			if(a.type == LEAF){
				//B becomes a leaf that is A but with the surface facing the opposite direction

				a.drawInfo.normal.multiplyScalar(-1);

				var i;
				for(i = 0; i < 12; i++){
					if(a.drawInfo.normals[i]){
						a.drawInfo.normals[i].multiplyScalar(-1);
					}
				}

				for(i = 0; i < 8; i++){
					if(a.drawInfo.corners[i] != MATERIAL_AIR){
						a.drawInfo.corners[i] = MATERIAL_AIR;
					}else{
						a.drawInfo.corners[i] = b.material;
					}
				}

				return a;
			}else if(a.type == INTERIOR){
				var i;
				for(i = 0; i < 8; i++){
					a.children[i] = this.getDifference(a.children[i],b);
				}
				return a;
			}
		}else{
			return b;
		}
	}


	if(a.type == LEAF && b.type == LEAF){
		var i;
		//If a corner in a is non air then that corner in b is air
		for(i = 0; i < 8; i++){
			if(a.drawInfo.corners[i] != MATERIAL_AIR){
				b.drawInfo.corners[i] = MATERIAL_AIR;
			}
		}


		var qefNeedsUpdate = false;
		//Edges in a override edges in b iff the intersection in a is closer to a non-air material
		for(i = 0; i < 12; i++){
			var c0 = edgevmap[i][0];
			var c1 = edgevmap[i][1];
			var edgeA = a.drawInfo.normals[i];
			var edgeB = b.drawInfo.normals[i];

			if(b.drawInfo.corners[c0] == MATERIAL_AIR && b.drawInfo.corners[c1] == MATERIAL_AIR){
					if(edgeB){
						//b.drawInfo.normal.sub(b.drawInfo.normals[i]);
						b.drawInfo.positions[i] = null;
						b.drawInfo.normals[i] = null; 
						qefNeedsUpdate = true;
					}
			}else if(!edgeB && edgeA){
				//B must be a solid edge of 1 material so just add intersection data of A with flipped normal
				b.drawInfo.normals[i] = a.drawInfo.normals[i];
				b.drawInfo.normals[i].multiplyScalar(-1);
				b.drawInfo.normal.add(b.drawInfo.normals[i]);
				b.drawInfo.positions[i] = a.drawInfo.positions[i];
				qefNeedsUpdate = true;

			}else if(edgeB && edgeA){
				//Find non air material corner in b.
				//If edge in a is closer than edge in B then override
				var solidPoint = new THREE.Vector3();
				if(b.drawInfo.corners[c0] != MATERIAL_AIR){
					solidPoint.addVectors(b.min,childOffsets[c0]*b.size);
				}else{
					solidPoint.addVectors(b.min,childOffsets[c1]*b.size);
				}

				if(a.drawInfo.positions[i].distanceTo(solidPoint) < b.drawInfo.positions[i].distanceTo(solidPoint)){
					b.drawInfo.positions[i] = a.drawInfo.positions[i];
					b.drawInfo.normals[i] = a.drawInfo.normals[i];
					b.drawInfo.normals[i].multiplyScalar(-1);
					b.drawInfo.normal.add(b.drawInfo.normals[i]);
					qefNeedsUpdate = true;
				}

			}
		}

		if(qefNeedsUpdate){
			b.drawInfo.qef.initialize(0,b.drawInfo.normals,b.drawInfo.positions);
			var qefPosition = new THREE.Vector3();
			b.drawInfo.qef.solve(qefPosition,QEF_ERROR, QEF_SWEEPS, QEF_ERROR);

			if(qefPosition.x < b.min.x || qefPosition.x > (b.min.x + b.size) ||
			qefPosition.y < b.min.y || qefPosition.y > (b.min.y + b.size) ||
			qefPosition.z < b.min.z || qefPosition.z > (b.min.z + b.size)){
				var mp = b.drawInfo.qef.getMassPoint();
				qefPosition = new THREE.Vector3(mp.x,mp.y,mp.z);
			}

			b.drawInfo.position = qefPosition;
		}

		b.drawInfo.normal.normalize();
		return b;
	}

}

extend(CSGNode,DiffNode);

//Intersect -> The volume of the left operand only where the right operand exists as well

IntersectNode = function(params){
	CSGNode.call(this,params);
}

extend (CSGNode,IntersectNode);

//Operations that produce a volume (set the world field and return self)
//Usually have a material index param

//Params:
//A world to set this node to
LoadNode = function(params){
	CSGNode.call(this,params);
}
LoadNode.prototype.operation = function(){
	var l = new CSGNode();
	l.world = this.params.world;
	return l;
}
extend(CSGNode,LoadNode);


//For an infinite world system we would have to run generation on each chunk that is affected
//For now we will say that a world is just one chunk

//Maybe these operations shouldn't be doing simplification?

//Params:
//THREE.Vector3 box origin
//THREE.Vector3 box halfDim
CuboidNode = function(params){
	CSGNode.call(this,params);
}

CuboidNode.prototype.operation = function(){
	var origin = this.params.origin;
	var halfDim = this.params.halfDim;

	var mat = 1;
	if(this.params.material){
		mat = this.params.material;
	}

	//How many samples of the density function should be taken
	var sampleSize = CHUNK_SIZE;
	if(this.params.sampleSize){
		sampleSize = this.params.sampleSize;
	}

	var cubeAABB = new AABB(new THREE.Vector3().subVectors(origin,halfDim),new THREE.Vector3(halfDim.x*2,halfDim.y*2,halfDim.z*2));
	var cubeFn = new cuboidDensityFn(sampleSize, cubeAABB, {origin: origin, halfDim: halfDim, material: mat});


	if(sampleSize <= CUBOID_CACHE_LIMIT){
		//Precompute density values
		cubeFn.cacheMaterial();
	}

	//For now we will only generate a chunk at the origin and if the operation is completely outside of it then it is invalid
	var c = new Chunk(new THREE.Vector3(0,0,0),cubeFn);

	c.buildOctree(c.root,false);
	//c.root = c.simplifyOctree(c.root,simplificationThreshold);

	var cNode = new CSGNode();
	cNode.world = c;
	return cNode;
}

extend(CSGNode,CuboidNode);

//Params:
//THREE.Vector3 sphere origin
//THREE.Vector3 sphere radius
//World (list of OctreeNodes representing)
SphereNode = function(params){
	CSGNode.call(this,params);
}

SphereNode.prototype.operation = function(){
	var origin = this.params.origin;
	var radius = this.params.radius;

	var mat = 1;
	if(this.params.material){
		mat = this.params.material;
	}


	//How many samples of the density function should be taken
	var sampleSize = CHUNK_SIZE;
	if(this.params.sampleSize){
		sampleSize = this.params.sampleSize;
	}


	var sphereAABB = new AABB(new THREE.Vector3().subVectors(origin,new THREE.Vector3(radius,radius,radius)),new THREE.Vector3(radius*2,radius*2,radius*2));
	var sphereFn = new sphereDensityFn(sampleSize, sphereAABB, {origin: origin, radius: radius, material: mat});

	if(sampleSize <= SPHERE_CACHE_LIMIT){
		//Precompute density values
		sphereFn.cacheMaterial();
	}

	//For now we will only generate a chunk at the origin and if the operation is completely outside of it then it does nothing
	var c = new Chunk(new THREE.Vector3(0,0,0),sphereFn);

	c.buildOctree(c.root,false);

	var sNode = new CSGNode();
	sNode.world = c;
	return sNode;
}
extend(CSGNode,SphereNode);
//Params
//THREE.Vector3 sample min
//THREE.Vector3 sample max
//World (list of OctreeNodes representing chunks of terrain)
PerlinHillsNode = function(params){
	CSGNode.call(this,params);
}
PerlinHillsNode.prototype.operation = function(){

	var grassDepth = 1;
	if(this.params.grassDepth){
		grassDepth = this.params.grassDepth;
	}

	var dirtDepth = 5;
	if(this.params.dirtDepth){
		dirtDepth = this.params.dirtDepth;
	}

	var sampleSize = CHUNK_SIZE;
	if(this.params.sampleSize){
		sampleSize = this.params.sampleSize;
	}

	var octaves = 4;
	if(this.params.octaves){
		octaves = this.params.octaves;
	}

	var frequency = 0.5;
	if(this.params.frequency){
		frequency = this.params.frequency;
	}

	var lacunarity = 2.2;
	if(this.params.lacunarity){
		lacunarity = this.params.lacunarity;
	}

	var persistence = 0.68;
	if(this.params.persistence){
		persistence = this.params.persistence;
	}

	var maxHeight = 8;
	if(this.params.maxHeight){
		maxHeight = this.params.maxHeight;
	}

	var scale = 64;
	if(this.params.scale){
		scale = this.params.scale;
	}

	var hillsFn = new perlinHillsDensityFn(sampleSize,{grassDepth: grassDepth,
														dirtDepth: dirtDepth,
														scale: scale,
														octaves: octaves,
														frequency: frequency,
														lacunarity: lacunarity,
														persistence: persistence,
														maxHeight: maxHeight});

	if(sampleSize <= NOISE_CACHE_LIMIT){
		hillsFn.cacheMaterial();
	}

	var c = new Chunk(new THREE.Vector3(0,0,0),hillsFn);

	c.buildOctree(c.root,true);

	var hNode = new CSGNode();
	hNode.world = c;
	return hNode;
}

extend(CSGNode,PerlinHillsNode);




PerlinCavesNode = function(params){
	CSGNode.call(this,params);
}

PerlinCavesNode.prototype.operation = function(){

	var sampleSize = CHUNK_SIZE;
	if(this.params.sampleSize){
		sampleSize = this.params.sampleSize;
	}

	var octaves = 4;
	if(this.params.octaves){
		octaves = this.params.octaves;
	}

	var frequency = 0.5;
	if(this.params.frequency){
		frequency = this.params.frequency;
	}

	var lacunarity = 2.2;
	if(this.params.lacunarity){
		lacunarity = this.params.lacunarity;
	}

	var persistence = 0.68;
	if(this.params.persistence){
		persistence = this.params.persistence;
	}

	var maxHeight = 8;
	if(this.params.maxHeight){
		maxHeight = this.params.maxHeight;
	}

	var scale = 64;
	if(this.params.scale){
		scale = this.params.scale;
	}

	var material = MATERIAL_DIRT;
	if(this.params.material){
		material = this.params.material;
	}

	var threshold = 0;
	if(this.params.threshold){
		threshold = this.params.threshold;
	}
	var cavesFn = new perlinCavesDensityFn(sampleSize,{material: material,
														scale: scale,
														octaves: octaves,
														frequency: frequency,
														lacunarity: lacunarity,
														persistence: persistence,
														threshold: threshold});

	if(sampleSize <= NOISE_CACHE_LIMIT){
		cavesFn.cacheMaterial();
	}

	var c = new Chunk(new THREE.Vector3(0,0,0),cavesFn);

	c.buildOctree(c.root,true);

	var hNode = new CSGNode();
	hNode.world = c;
	return hNode;
}

extend(CSGNode,PerlinCavesNode);



densityFunction = function(sampleSize, AABB, params){
	this.sampleSize = sampleSize;

	this.leafSize = CHUNK_SIZE / sampleSize;
	this.indexConversion = 1/this.leafSize;

	this.AABB = AABB;
	this.params = params;
};


densityFunction.prototype.cacheMaterial = function(){
	this.material = [];
	var x,y,z;

	for(x = 0; x < this.sampleSize+1; x++){
		this.material[x] = [];
		for(y = 0; y < this.sampleSize+1; y++){
			this.material[x][y] = [];
			for(z = 0; z < this.sampleSize+1; z++){
				this.getDensity(x*this.leafSize,y*this.leafSize,z*this.leafSize);
				this.material[x][y][z] = this.dominantMaterial;
			}
		}
	}
};

densityFunction.prototype.getMaterialData = function(x,y,z){
	if(this.material){
		return this.material[x*this.indexConversion][y*this.indexConversion][z*this.indexConversion];
	}else{
		this.getDensity(x,y,z);
		return this.dominantMaterial;
	}
}

//Calculate the normal of the density function using the finite difference method
//Could probably get a better normal by calculating the gradient analytically... but it wouldn't be faster
densityFunction.prototype.getSurfaceNormal = function(v){
	var x = v.x;
	var y = v.y;
	var z = v.z;
	var h = 0.001;
	var dx = this.getDensity(x+h,y,z) - this.getDensity(x-h,y,z);
	var dy = this.getDensity(x,y+h,z) - this.getDensity(x,y-h,z);
	var dz = this.getDensity(x,y,z+h) - this.getDensity(x,y,z-h);
	var normal = new THREE.Vector3(dx,dy,dz).normalize();
	return normal;
};

densityFunction.prototype.getSurfaceIntersection = function(p1, p2){
	var minVal = 100000;
	var t = 0;
	var currentT = 0;
	var steps = 8.0;
	var increment = 1.0/steps;

	while(currentT <= 1.0){
		var p = new THREE.Vector3();
		var l = new THREE.Vector3();
		l.subVectors(p2,p1).multiplyScalar(currentT);
		p.addVectors(p1,l);
		var density = Math.abs(this.getDensity(p.x,p.y,p.z));
		if(density < minVal){
			minVal = density;
			t = currentT;
		}

		currentT += increment;
	}

	l.subVectors(p2,p1).multiplyScalar(t);
	var result = new THREE.Vector3();
	result.addVectors(p1,l);
	return result;
};

//Different kinds of density functions
cuboidDensityFn = function(sampleSize, AABB, params){
	densityFunction.call(this,sampleSize,AABB,params);
}

cuboidDensityFn.prototype.getDensity = function(x,y,z){
	var worldPosition = new THREE.Vector3(x,y,z);
	var localPosition = new THREE.Vector3().subVectors(worldPosition,this.params.origin);
	var abso = new THREE.Vector3(Math.abs(localPosition.x),Math.abs(localPosition.y),Math.abs(localPosition.z));
	var d = new THREE.Vector3().subVectors(abso,this.params.halfDim);
	var m = Math.max(d.x,Math.max(d.y,d.z));

	var density = Math.min(m, d.max(new THREE.Vector3(0,0,0)).length());

	if(density < 0){
		this.dominantMaterial = this.params.material;
	}else{
		this.dominantMaterial = MATERIAL_AIR;
	}

	return density;
	
};

extend(densityFunction,cuboidDensityFn);



sphereDensityFn = function(sampleSize, AABB, params){
	densityFunction.call(this,sampleSize,AABB,params);
}

sphereDensityFn.prototype.getDensity = function(x,y,z){
	var worldPosition = new THREE.Vector3(x,y,z);


	var density = worldPosition.sub(this.params.origin).length() - this.params.radius;

	if(density < 0){
		this.dominantMaterial = this.params.material;
	}else{
		this.dominantMaterial = MATERIAL_AIR;
	}

	return density
};

extend(densityFunction,sphereDensityFn);



perlinHillsDensityFn = function(sampleSize, params){
	densityFunction.call(this,sampleSize,null,params);
}

//Only need to cache a heightmap
perlinHillsDensityFn.prototype.cacheMaterial = function(){
	this.material = [];
	var x,y;
	for(x = 0; x < this.sampleSize+1; x++){
		this.material[x] = [];
		for(y = 0; y < this.sampleSize+1; y++){
			this.material[x][y] = this.getHeight(x,y);
		}
	}
}

perlinHillsDensityFn.prototype.getSurfaceIntersection = function(p1,p2){
	var minVal = 100000;
	var t = 0;
	var steps = 8.0;
	var increment = 1.0/steps;
	var currentT = increment;
	var densityOffset = 0;

	var p1Density = this.getDensity(p1.x,p1.y,p1.z);
	var p2Density = this.getDensity(p2.x,p2.y,p2.z);


	if(p1Density < -this.dirtDepth || p2Density < -this.dirtDepth){
		densityOffset = this.dirtDepth;
	}else if(p1Density < -this.grassDepth || p2Density < -this.grassDepth){
		densityOffset = this.grassDepth;
	}

	if(Math.abs(p1Density+densityOffset) < Math.abs(p2Density+densityOffset)){
		minVal = Math.abs(p1Density+densityOffset);
		t = 0;
	}else{
		minVal = Math.abs(p2Density+densityOffset);
		t = 1;
	}

	while(currentT <= 1.0-increment){
		var p = new THREE.Vector3();
		var l = new THREE.Vector3();
		l.subVectors(p2,p1).multiplyScalar(currentT);
		p.addVectors(p1,l);
		var density = Math.abs(this.getDensity(p.x,p.y,p.z)+densityOffset);
		if(density < minVal){
			minVal = density;
			t = currentT;
		}

		currentT += increment;
	}

	l.subVectors(p2,p1).multiplyScalar(t);
	var result = new THREE.Vector3();
	result.addVectors(p1,l);
	return result;
}

perlinHillsDensityFn.prototype.getSurfaceNormal = function(v){
	var samplePoint = new THREE.Vector3(v.x,this.getHeight(v.x,v.z),v.z);
	return densityFunction.prototype.getSurfaceNormal.call(this,samplePoint);
}

perlinHillsDensityFn.prototype.getDensity = function(x,y,z){
	var  height = this.getHeight(x,z);
	return y - height;
}

perlinHillsDensityFn.prototype.getHeight = function(x,y){
	var SCALE = 1 / this.params.scale;

	x*=SCALE; 
	y*=SCALE;

	var height = 0;

	var amplitude = 1;
	x*=this.params.frequency;
	y*=this.params.frequency;

	var i;
	for(i = 0; i < this.params.octaves; i++){
		height += noise.perlin2(x,y)*amplitude;
		x *= this.params.lacunarity;
		y *= this.params.lacunarity;
		amplitude *= this.params.persistence;
	}

	height = (0.5+height*0.5)*this.params.maxHeight;
	if(height < 0){
		height = 0.1;
	}
	return height;
}

perlinHillsDensityFn.prototype.getMaterialData = function(x,y,z){

	var density;
	if(this.material){
		density = y - this.material[x*this.indexConversion][z*this.indexConversion];
	}else{
		density = this.getDensity(x,y,z);
	}
	

	if(density < -this.params.dirtDepth){
		return MATERIAL_ROCK;
	}else if(density < -this.params.grassDepth){
		return MATERIAL_DIRT;
	}else if(density < 0){
		return MATERIAL_GRASS;
	}else{
		return MATERIAL_AIR;
	}
}

extend(densityFunction, perlinHillsDensityFn);





perlinCavesDensityFn = function(sampleSize, params){
	densityFunction.call(this,sampleSize,null,params);
	this.densityThreshold = params.threshold;
}

perlinCavesDensityFn.prototype.getDensity = function(x,y,z){
	var density = 0;

	var SCALE = 1 / this.params.scale;

	x*=SCALE; 
	y*=SCALE;
	z*=SCALE;

	var amplitude = 1;
	x*=this.params.frequency;
	y*=this.params.frequency;
	z*=this.params.frequency;

	var normalization = 0;

	var i;
	for(i = 0; i < this.params.octaves; i++){
		normalization += amplitude;
		density += noise.simplex3(x,z,y)*amplitude;
		x *= this.params.lacunarity;
		y *= this.params.lacunarity;
		z *= this.params.lacunarity;
		amplitude *= this.params.persistence;
	}
	density -= this.densityThreshold;

	density;


	if(density > this.params.threshold){
		this.dominantMaterial = this.params.material;
	}else{
		this.dominantMaterial = MATERIAL_AIR;
	}

	return density;
}

perlinCavesDensityFn.prototype.getSurfaceNormal = function(v){
	return densityFunction.prototype.getSurfaceNormal.call(this,v).multiplyScalar(-1);
}

extend(densityFunction, perlinCavesDensityFn);








AABB = function(min,dim){
	this.min = min;
	this.max = new THREE.Vector3().addVectors(min,dim);
};

AABB.prototype.intersects = function(other){
	if(this.max.x < other.min.x)
		return false;
	if(this.min.x > other.max.x)
		return false;
	if(this.max.y < other.min.y)
		return false;
	if(this.min.y > other.max.y)
		return false;
	if(this.max.z < other.min.z)
		return false;
	if(this.min.z > other.max.z)
		return false;
	return true;
};

AABB.prototype.contains = function(other){
	if(other.min.x < this.min.x)
		return false;
	if(other.max.x > this.max.x)
		return false;
	if(other.min.y < this.min.y)
		return false;
	if(other.max.y > this.max.y)
		return false;
	if(other.min.z < this.min.z)
		return false;
	if(other.max.z > this.max.z)
		return false;
	return true;
};




