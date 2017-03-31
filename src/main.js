
var scene, camera, renderer, controls;
var geometry, material, mesh;
var clock = new THREE.Clock();



document.documentElement.style.overflow = 'hidden';
document.body.scroll = "no";
document.getElementsByTagName("body")[0].style.marginLeft = 0;
document.getElementsByTagName("body")[0].style.marginTop = 0;

var terrainMesh;
var terrainNode;
var cubeBrushDim = new THREE.Vector3(1,1,1);
var sphereBrushRadius = 1;

var brushShape = "sphere";

var brushSizeDelta = 0.25;

var brushLocation = new THREE.Vector3(32,40,32);
var brushPositionDelta = 0.25;

var brushMaterial = MATERIAL_AIR;

var sphereBrush;
var cubeBrush;



init();
animate();
simulate();

/* Get a THREE.js Geometry object from the root node of a world */
function getTerrainGeometry(terrainNode, simplificationThreshold){

	var terrain = terrainNode.world;
	terrain.generateSimplifiedMesh(simplificationThreshold);

	var terrainGeometry = new THREE.Geometry();
	var i;

	for(i = 0; i < terrain.vertexCount; i++){

		var vPos = terrain.vertexBuffer[i].xyz;
		terrainGeometry.vertices.push(vPos);
	}

	var faceCount = 0;
	for(i = 0; i < terrain.indexCount; i += 3){
		var a = terrain.indexBuffer[i];
		var b = terrain.indexBuffer[i+1];
		var c = terrain.indexBuffer[i+2];

		var face = new THREE.Face3(a,b,c);

		face.vertexNormals[0] = terrain.vertexBuffer[a].normal;
		face.vertexNormals[1] = terrain.vertexBuffer[b].normal;
		face.vertexNormals[2] = terrain.vertexBuffer[c].normal;
		face.normal = new THREE.Vector3().addVectors(terrain.vertexBuffer[terrain.indexBuffer[i]].normal,
			new THREE.Vector3().addVectors(terrain.vertexBuffer[terrain.indexBuffer[i+1]].normal,
				terrain.vertexBuffer[terrain.indexBuffer[i+2]].normal));
		face.normal.normalize();

		face.materialIndex = terrain.materialBuffer[faceCount];

		terrainGeometry.faces.push(face);
		faceCount++;
	}

	return terrainGeometry;
}

function init(){
	scene = new THREE.Scene();


	camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 1, 10000 );


	renderer = new THREE.WebGLRenderer();
	renderer.setSize(window.innerWidth, window.innerHeight);

	window.addEventListener('resize',onWindowResize, false);
	document.body.appendChild(renderer.domElement);


	var fov = 80;

	controls = new THREE.TrackballControls(camera,renderer.domElement);
	controls.target.set(32,16,32);

	camera.position.y = 100;

	noise.seed(Math.random());

	var origin = new THREE.Vector3(0,0,0);


	var wall1 = new CuboidNode({material: MATERIAL_ROCK, 
								origin: new THREE.Vector3(0,32,32), 
								halfDim: new THREE.Vector3(1,32,32),
								sampleSize: CHUNK_SIZE});


	var wall2 = new CuboidNode({material: MATERIAL_ROCK, 
								origin: new THREE.Vector3(32,32,0), 
								halfDim: new THREE.Vector3(32,32,1),
								sampleSize: CHUNK_SIZE});


	var wall3 = new CuboidNode({material: MATERIAL_ROCK, 
								origin: new THREE.Vector3(64,32,32), 
								halfDim: new THREE.Vector3(1,32,32),
								sampleSize: CHUNK_SIZE});

	var wall4 = new CuboidNode({material: MATERIAL_GRASS, 
								origin: new THREE.Vector3(32,32,64), 
								halfDim: new THREE.Vector3(32,32,1),
								sampleSize: CHUNK_SIZE});

	var floor = new CuboidNode({material: MATERIAL_GRASS, 
								origin: new THREE.Vector3(32,0,32), 
								halfDim: new THREE.Vector3(32,1,32),
								sampleSize: CHUNK_SIZE});


	var hillsNode = new PerlinHillsNode({grassDepth: 2,
										dirtDepth: 4,
										maxHeight: 50,
										octaves: 4,
										frequency: 0.3,
										lacunarity: 2.2,
										persistence: 0.68,
										sampleSize: CHUNK_SIZE});

	var caveNode = new PerlinCavesNode({material: MATERIAL_ROCK,
										octaves: 2,
										frequency: 0.65,
										lacunarity: 1, //Connectedness
										persistence: 0.3,
										scale: 16,
										threshold: 0.3,
										sampleSize: CHUNK_SIZE});



	var unionNode = new UnionNode();
	unionNode.addChild(wall1);
	unionNode.addChild(wall2);
	unionNode.addChild(wall3);
	unionNode.addChild(wall4);
	unionNode.addChild(floor);



	var diffNode = new DiffNode();
	diffNode.addChild(caveNode);
	diffNode.addChild(hillsNode);

	terrainNode = new DiffNode();
	terrainNode.addChild(unionNode);
	terrainNode.addChild(diffNode);

	terrainNode = terrainNode.getResult();

	var materials = [];
	materials[0] = new THREE.MeshLambertMaterial({color: 0xff0000});
	materials[MATERIAL_DIRT] = new THREE.MeshLambertMaterial({color: 0x6d3f17});
	materials[MATERIAL_ROCK] = new THREE.MeshLambertMaterial({color: 0x7e838c});
	materials[MATERIAL_GRASS] = new THREE.MeshLambertMaterial({color: 0x186813});

	var material = new THREE.MultiMaterial(materials);

	terrainMesh = new THREE.Mesh(getTerrainGeometry(terrainNode,0.3), material);


	unionNode = null;
	diffNode = null;
	wall1 = null;
	wall2 = null;
	wall3 = null;
	wall4 = null;
	floor = null;
	hillsNode = null;
	caveNode = null;

	//terrain.drawOctree(terrain.root,scene);

	scene.add(terrainMesh);


	var directionalLight = new THREE.DirectionalLight( 0xffffff, 1 );
	directionalLight.position.set( 0, 1, 0 );
	scene.add( directionalLight );

	var light = new THREE.AmbientLight( 0x404040 );
	scene.add( light );

	var brushMaterial = new THREE.MeshBasicMaterial({color: 0xdddd00, transparent: true, opacity: 0.7});

	sphereBrush = new THREE.Mesh(new THREE.SphereGeometry(1,32,32),brushMaterial);
	scene.add(sphereBrush);
	sphereBrush.position.set(brushLocation.x,brushLocation.y,brushLocation.z);

	cubeBrush = new THREE.Mesh(new THREE.BoxGeometry(1,1,1),brushMaterial);
	scene.add(cubeBrush);
	cubeBrush.visible = false;

	document.addEventListener("keydown",onKeyDown, false);
	document.addEventListener("keyup",onKeyUp,false);
}

var useWireframe = false;
function onKeyDown(event){
	var keyCode = event.which;

	if(keyCode == 79){//Radius/Y increase
		if(brushShape === "sphere"){
			sphereBrushRadius += brushSizeDelta;
			sphereBrush.scale.set(sphereBrushRadius,sphereBrushRadius,sphereBrushRadius);
		}else if(brushShape === "box"){
			cubeBrushDim.set(cubeBrushDim.x,cubeBrushDim.y+brushSizeDelta,cubeBrushDim.z);
			cubeBrush.scale.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z);
		}

	}else if(keyCode == 76){//Radius/Y decrease
		if(brushShape === "sphere"){
			sphereBrushRadius -= brushSizeDelta;
			if(sphereBrushRadius < 1){
				sphereBrushRadius = 1;
			}
			sphereBrush.scale.set(sphereBrushRadius,sphereBrushRadius,sphereBrushRadius);
		}else if(brushShape === "box"){
			cubeBrushDim.set(cubeBrushDim.x,cubeBrushDim.y-brushSizeDelta,cubeBrushDim.z);
			if(cubeBrushDim.y < 1){
				cubeBrushDim.set(cubeBrushDim.x,1,cubeBrushDim.z);
			}
			cubeBrush.scale.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z);
		}
	}else if(keyCode == 90){//wireframe toggle
		useWireframe = !useWireframe;
		terrainMesh.material.materials[1].wireframe = useWireframe;
		terrainMesh.material.materials[2].wireframe = useWireframe;
		terrainMesh.material.materials[3].wireframe = useWireframe;
	}else if(keyCode == 84){//Toggle brush shape
		if(brushShape === "sphere"){
			brushShape = "box";
			cubeBrush.visible = true;
			sphereBrush.visible = false;
			cubeBrushDim.set(sphereBrushRadius,sphereBrushRadius,sphereBrushRadius);
			cubeBrush.scale.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z);
			cubeBrush.position.set(brushLocation.x,brushLocation.y,brushLocation.z);
		}else{
			brushShape = "sphere";
			cubeBrush.visible = false;
			sphereBrush.visible = true;
			sphereBrush.position.set(brushLocation.x,brushLocation.y,brushLocation.z);
		}
	}else if(keyCode == 73){//z+
			cubeBrushDim.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z+brushSizeDelta);
			cubeBrush.scale.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z);
	}else if(keyCode == 77){//z-
			cubeBrushDim.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z-brushSizeDelta);
			if(cubeBrushDim.z < 1){
				cubeBrushDim.set(cubeBrushDim.x,cubeBrushDim.y,1);
			}
			cubeBrush.scale.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z);
	}else if(keyCode == 74){//x-
			cubeBrushDim.set(cubeBrushDim.x-brushSizeDelta,cubeBrushDim.y,cubeBrushDim.z);
			if(cubeBrushDim.x < 1){
				cubeBrushDim.set(1,cubeBrushDim.y,cubeBrushDim.z);
			}
			cubeBrush.scale.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z);
	}else if(keyCode == 75){//x+
			cubeBrushDim.set(cubeBrushDim.x+brushSizeDelta,cubeBrushDim.y,cubeBrushDim.z);
			cubeBrush.scale.set(cubeBrushDim.x,cubeBrushDim.y,cubeBrushDim.z);
	}else if(keyCode == 89){//pos z+
		brushLocation.set(brushLocation.x,brushLocation.y,brushLocation.z+brushPositionDelta);
	}else if(keyCode == 66){//pos z-
		brushLocation.set(brushLocation.x,brushLocation.y,brushLocation.z-brushPositionDelta);
	}else if(keyCode == 72){//pos x+
		brushLocation.set(brushLocation.x+brushPositionDelta,brushLocation.y,brushLocation.z);
	}else if(keyCode == 71){//pos x-
		brushLocation.set(brushLocation.x-brushPositionDelta,brushLocation.y,brushLocation.z);
	}else if(keyCode == 85){//pos y+
		brushLocation.set(brushLocation.x,brushLocation.y+brushPositionDelta,brushLocation.z);
	}else if(keyCode == 78){//pos y-
		brushLocation.set(brushLocation.x,brushLocation.y-brushPositionDelta,brushLocation.z);
	}else if(keyCode == 49){
		brushMaterial = MATERIAL_AIR;
	}else if(keyCode == 50){
		brushMaterial = MATERIAL_DIRT;
	}else if(keyCode == 51){
		brushMaterial = MATERIAL_ROCK
	}else if(keyCode == 52){
		brushMaterial = MATERIAL_GRASS;
	}

	cubeBrush.position.set(brushLocation.x,brushLocation.y,brushLocation.z);
	sphereBrush.position.set(brushLocation.x,brushLocation.y,brushLocation.z);
}

function onKeyUp(event){
	var keyCode = event.which;

	if(keyCode == 13){
		if(brushMaterial == MATERIAL_AIR){
			var diffNode = new DiffNode();
			var subtractNode;

			if(brushShape === "sphere"){
				subtractNode = new SphereNode({material: MATERIAL_DIRT,
												origin: brushLocation,
												radius: sphereBrushRadius,
												sampleSize: CHUNK_SIZE});
			}else{
				subtractNode = new CuboidNode({material: MATERIAL_DIRT,
												origin: brushLocation,
												halfDim: new THREE.Vector3(cubeBrush.scale.x/2,cubeBrush.scale.y/2,cubeBrush.scale.z/2),
												sampleSize: CHUNK_SIZE});
			}


			diffNode.addChild(subtractNode);
			diffNode.addChild(terrainNode);

			terrainNode = diffNode.getResult();

			terrainMesh.geometry.dispose();
			terrainMesh.geometry = getTerrainGeometry(terrainNode,0.3).clone();


		}else{
			var unionNode = new UnionNode();
			var addNode;

			if(brushShape === "sphere"){
				addNode = new SphereNode({material: brushMaterial,
												origin: brushLocation,
												radius: sphereBrushRadius,
												sampleSize: CHUNK_SIZE});
			}else{
				addNode = new CuboidNode({material: brushMaterial,
												origin: brushLocation,
												halfDim: new THREE.Vector3(cubeBrush.scale.x/2,cubeBrush.scale.y/2,cubeBrush.scale.z/2),
												sampleSize: CHUNK_SIZE});
			}

			unionNode.addChild(terrainNode);
			unionNode.addChild(addNode);

			terrainNode = unionNode.getResult();

			terrainMesh.geometry.dispose();
			terrainMesh.geometry = getTerrainGeometry(terrainNode,0.3).clone();
		}
	}
}

function onWindowResize(){
	camera.aspect = window.innerWidth / window.innerHeight;
	camera.updateProjectionMatrix();

	renderer.setSize(window.innerWidth, window.innerHeight);
}

function animate(){
	var delta = clock.getDelta();
	requestAnimationFrame(animate);
	controls.update(delta);
	renderer.render(scene, camera);
}

function simulate(){

	setTimeout(simulate, 1000 / 35);
}