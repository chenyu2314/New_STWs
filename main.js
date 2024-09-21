
// This code is based on three.js, which comes with the following license:
//
// The MIT License
//
// Copyright © 2010-2024 three.js authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
import * as THREE from 'three';

import { GUI } from 'three/addons/libs/lil-gui.module.min.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

let name = 'LensletLegend';
let scene;
let aspectRatioVideoFeedU = 4.0/3.0;
let aspectRatioVideoFeedE = 4.0/3.0;
let renderer;
let videoFeedU, videoFeedE;	// feeds from user/environment-facing cameras
let camera;
let controls;
let raytracingSphere;
let raytracingSphereShaderMaterial;
	
// Nokia HR20, according to https://www.camerafv5.com/devices/manufacturers/hmd_global/nokia_xr20_ttg_0/
let fovVideoFeedU = 67.3;	// (user-facing) camera
let fovVideoFeedE = 68.3;	// (environment-facing) camera
let fovScreen = 68;

let cameraLensDistance = 10.0;
let raytracingSphereRadius = 20.0;
let offsetFromConfocal = 0.0;
let deltaPeriod = 0.0;

// camera with wide aperture
let apertureRadius = 0.0;
let focusDistance = 1e8;
let noOfRays = 1;
// let pointsOnAperture = [];


// the status text area
let status = document.createElement('div');
let statusTime;	// the time the last status was posted

// the info text area
let info = document.createElement('div');
let infoTime;	// the time the last info was posted

let gui;

// let counter = 0;

// true if stored photo is showing
let showingStoredPhoto = false;
let storedPhoto;
let storedPhotoDescription;
let storedPhotoInfoString;


//
let fMax = 0.5;
let magnification = 0.5;
let rotationAngle = 0;
let gapBetween2And3 = fMax*magnification; let gapBetween4And5 = fMax*magnification*magnification;
let distance12 = 0.0001;let distance34 = 0.0001;let distance56 = 0.0001;
let basePeriod = 0.4;
//


// my Canon EOS450D
const click = new Audio('./click.m4a');

// uncomment the foolowing lines, and 
// stats.begin();
// stats.end();
// in animate(), to show fps stats
// import Stats from 'stats.js'
// var stats = new Stats();
// stats.showPanel( 0 ); // 0: fps, 1: ms, 2: mb, 3+: custom
// document.body.appendChild( stats.dom );

init();
animate();

function init() {
	// create the info element first so that any problems can be communicated
	createStatus();

	scene = new THREE.Scene();
	// scene.background = new THREE.Color( 'skyblue' );
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera = new THREE.PerspectiveCamera( fovScreen, windowAspectRatio, 0.1, 2*raytracingSphereRadius + 1 );
	camera.position.z = cameraLensDistance;
	screenChanged();
	
	renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
	renderer.setPixelRatio(window.devicePixelRatio);
	renderer.setSize( window.innerWidth, window.innerHeight );
	document.body.appendChild( renderer.domElement );
	// document.getElementById('livePhoto').appendChild( renderer.domElement );

	createVideoFeeds();

	addRaytracingSphere();

	// user interface

	addEventListenersEtc();

	addOrbitControls();

	// the controls menu
	createGUI();

	createInfo();
	refreshInfo();
}

function animate() {
	requestAnimationFrame( animate );

	// stats.begin();

	if(!showingStoredPhoto) {
		// update uniforms
		updateUniforms();

		renderer.render( scene,  camera );
	}

	// stats.end();
}

function updateUniforms() {


	// 控制透镜对之间的距离
    let distances = [distance12, distance34, distance56];  // 控制1号-2号、3号-4号、5号-6号之间的距离
    let gaps = [gapBetween2And3, gapBetween4And5];  // 控制2号-3号、4号-5号之间的间隔

    // 累加变量，控制透镜沿 z 轴的位置
    let currentZ = -1;
	
    // 循环更新每个透镜的中心位置、旋转角度和其他属性
    for (let i = 0; i < 6; i++) {
        // 更新透镜的中心位置
        raytracingSphereShaderMaterial.uniforms.centreOfArray.value[i].z = currentZ;
		
        // 根据透镜的位置关系调整 z 坐标
        currentZ -= distances[Math.floor(i / 2)];

        // 处理2号与3号、4号与5号之间的间隔
        if (i === 1 || i === 3) {
            currentZ -= gaps[Math.floor(i / 2)];
        }

        // 更新透镜的旋转角度
        raytracingSphereShaderMaterial.uniforms.alpha.value[i] = rotationAngle;

        // 更新透镜的周期（period）
        if (i === 0) {
            raytracingSphereShaderMaterial.uniforms.period.value[i] = basePeriod;
        } else {
            raytracingSphereShaderMaterial.uniforms.period.value[i] = basePeriod + i * deltaPeriod;
        }

    }
	

	// the tangents for the environment-facing camera video feed
	let tanHalfFovHE, tanHalfFovVE;
	if(aspectRatioVideoFeedE > 1.0) {
		// horizontal orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)/aspectRatioVideoFeedE;
	} else {
		// vertical orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)*aspectRatioVideoFeedE;
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.halfWidthE.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovHE;
	raytracingSphereShaderMaterial.uniforms.halfHeightE.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovVE;

	// the tangents for the user-facing camera video feed
	let tanHalfFovHU, tanHalfFovVU;
	if(aspectRatioVideoFeedU > 1.0) {
		// horizontal orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)/aspectRatioVideoFeedU;
	} else {
		// vertical orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)*aspectRatioVideoFeedU;
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.halfWidthU.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovHU;
	raytracingSphereShaderMaterial.uniforms.halfHeightU.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovVU;

	// // calculate the separation between the two arrays, s = f1 + f2 + offsetFromConfocal
	// let s = 
	// 	raytracingSphereShaderMaterial.uniforms.lensletsF1.value + 
	// 	raytracingSphereShaderMaterial.uniforms.lensletsF2.value +
	// 	offsetFromConfocal;
	// // arrange them symmetrically around z=0
	// raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.z = 0;	// +0.5*s;
	// raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.z = -s;	// -0.5*s;	// - 0.0001;

	// set the array periods
	//raytracingSphereShaderMaterial.uniforms.period2.value = raytracingSphereShaderMaterial.uniforms.period1.value + deltaPeriod;

	// create the points on the aperture

	// create basis vectors for the camera's clear aperture
	let viewDirection = new THREE.Vector3();
	let apertureBasisVector1 = new THREE.Vector3();
	let apertureBasisVector2 = new THREE.Vector3();
	camera.getWorldDirection(viewDirection);
	// if(counter < 10) console.log(`viewDirection = (${viewDirection.x.toPrecision(2)}, ${viewDirection.y.toPrecision(2)}, ${viewDirection.z.toPrecision(2)})`);

	if((viewDirection.x == 0.0) && (viewDirection.y == 0.0)) {
		// viewDirection is along z direction
		apertureBasisVector1.crossVectors(viewDirection, new THREE.Vector3(1, 0, 0)).normalize();
	} else {
		// viewDirection is not along z direction
		apertureBasisVector1.crossVectors(viewDirection, new THREE.Vector3(0, 0, 1)).normalize();
	}
	// viewDirection = new THREE.Vector3(0, 0, -1);
	// apertureBasisVector1 = new THREE.Vector3(1, 0, 0);
	apertureBasisVector2.crossVectors(viewDirection, apertureBasisVector1).normalize();

	// apertureBasis1 *= apertureRadius;
	// apertureBasis2 *= apertureRadius;

	// if(counter < 10) console.log(`apertureBasisVector1 = (${apertureBasisVector1.x.toPrecision(2)}, ${apertureBasisVector1.y.toPrecision(2)}, ${apertureBasisVector1.z.toPrecision(2)})`);
	// if(counter < 10) console.log(`apertureBasisVector2 = (${apertureBasisVector2.x.toPrecision(2)}, ${apertureBasisVector2.y.toPrecision(2)}, ${apertureBasisVector2.z.toPrecision(2)})`);
	// counter++;

	// create random points on the (circular) aperture
	// let i=0;
	// pointsOnAperture = [];	// clear the array containing points on the aperture
	// do {
	// 	// create a new random point on the camera's clear aperture
	// 	let x = 2*Math.random()-1;	// random number between -1 and 1
	// 	let y = 2*Math.random()-1;	// random number between -1 and 1
	// 	if(x*x + y*y <= 1) {
	// 		// (x,y) lies within a circle of radius 1
	// 		//  add a new point to the array of points on the aperture
	// 		pointsOnAperture.push(apertureRadius*x*apertureBasis1 + apertureRadius*y*apertureBasis2);
	// 		i++;
	// 	}
	// } while (i < noOfRays);
	raytracingSphereShaderMaterial.uniforms.noOfRays.value = noOfRays;
	raytracingSphereShaderMaterial.uniforms.apertureXHat.value.x = apertureRadius*apertureBasisVector1.x;
	raytracingSphereShaderMaterial.uniforms.apertureXHat.value.y = apertureRadius*apertureBasisVector1.y;
	raytracingSphereShaderMaterial.uniforms.apertureXHat.value.z = apertureRadius*apertureBasisVector1.z;
	raytracingSphereShaderMaterial.uniforms.apertureYHat.value.x = apertureRadius*apertureBasisVector2.x;
	raytracingSphereShaderMaterial.uniforms.apertureYHat.value.y = apertureRadius*apertureBasisVector2.y;
	raytracingSphereShaderMaterial.uniforms.apertureYHat.value.z = apertureRadius*apertureBasisVector2.z;
	// raytracingSphereShaderMaterial.uniforms.pointsOnAperture.value = pointsOnAperture;
	raytracingSphereShaderMaterial.uniforms.apertureRadius.value = apertureRadius;
	raytracingSphereShaderMaterial.uniforms.focusDistance.value = focusDistance;
	
	// (re)create random numbers
	// let i=0;
	// let randomNumbersX = [];
	// let randomNumbersY = [];
	// do {
	// 	// create a new pairs or random numbers (x, y) such that x^2 + y^2 <= 1
	// 	let x = 2*Math.random()-1;	// random number between -1 and 1
	// 	let y = 2*Math.random()-1;	// random number between -1 and 1
	// 	if(x*x + y*y <= 1) {
	// 		// (x,y) lies within a circle of radius 1
	// 		//  add a new point to the array of points on the aperture
	// 		randomNumbersX.push(apertureRadius*x);
	// 		randomNumbersY.push(apertureRadius*y);
	// 		i++;
	// 	}
	// } while (i < 100);
	// raytracingSphereShaderMaterial.uniforms.randomNumbersX.value = randomNumbersX;
	// raytracingSphereShaderMaterial.uniforms.randomNumbersY.value = randomNumbersY;
}

function createVideoFeeds() {
	// create the video stream for the user-facing camera first, as some devices (such as my iPad), which have both cameras,
	// but can (for whatever reason) only have a video feed from one at a time, seem to go with the video stream that was
	// created last, and as the standard view is looking "forward" it is preferable to see the environment-facing camera.
	videoFeedU = document.getElementById( 'videoFeedU' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// user-facing camera
		const constraintsU = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'user'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsU ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedU.srcObject = stream;
			videoFeedU.play();

			videoFeedU.addEventListener("playing", () => {
				aspectRatioVideoFeedU = videoFeedU.videoWidth / videoFeedU.videoHeight;
				updateUniforms();
				postStatus(`User-facing(?) camera resolution ${videoFeedU.videoWidth} &times; ${videoFeedU.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access user-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}

	videoFeedE = document.getElementById( 'videoFeedE' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// environment-facing camera
		const constraintsE = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'environment'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsE ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedE.srcObject = stream;
			videoFeedE.play();

			videoFeedE.addEventListener("playing", () => {
				aspectRatioVideoFeedE = videoFeedE.videoWidth / videoFeedE.videoHeight;
				updateUniforms();
				postStatus(`Environment-facing(?) camera resolution ${videoFeedE.videoWidth} &times; ${videoFeedE.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access environment-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}
}

/** create raytracing phere */
function addRaytracingSphere() {
	const videoFeedUTexture = new THREE.VideoTexture( videoFeedU );
	const videoFeedETexture = new THREE.VideoTexture( videoFeedE );
	videoFeedUTexture.colorSpace = THREE.SRGBColorSpace;
	videoFeedETexture.colorSpace = THREE.SRGBColorSpace;

	// create arrays of random numbers (as GLSL is rubbish at doing random numbers)
	let randomNumbersX = [];
	let randomNumbersY = [];
	// make the first random number 0 in both arrays, meaning the 0th ray starts from the centre of the aperture
	randomNumbersX.push(0);
	randomNumbersY.push(0);
	// fill in the rest of the array with random numbers
	let i=1;
	do {
		// create a new pairs or random numbers (x, y) such that x^2 + y^2 <= 1
		let x = 2*Math.random()-1;	// random number between -1 and 1
		let y = 2*Math.random()-1;	// random number between -1 and 1
		if(x*x + y*y <= 1) {
			// (x,y) lies within a circle of radius 1
			//  add a new point to the array of points on the aperture
			randomNumbersX.push(x);
			randomNumbersY.push(y);
			i++;
		}
	} while (i < 100);

	// the sphere surrouning the camera in all directions
	const geometry = 
		new THREE.SphereGeometry( raytracingSphereRadius );
	raytracingSphereShaderMaterial = new THREE.ShaderMaterial({
		side: THREE.DoubleSide,
		// wireframe: true,
		uniforms: {
			displayMode: { value:0 },
			visible: { value: Array(6).fill(true) }, // visible array with 6 elements, all true
			period: { value: Array(6).fill(0.4) }, // period array with 6 elements, all 0.4
			alpha: { value: Array(6).fill(0.0) }, // alpha array with 6 elements, all 0.0
			lensletsF: { value: [fMax, -fMax, magnification*fMax, -magnification*fMax, 0, -magnification*magnification*fMax] }, // lensletsF array with 6 elements, all 1.0, Array(6).fill(1.0)
			additionalF: { value: Array(6).fill(1e10) }, // additionalF array with 6 elements, all 1e10
			centreOfArray: { 
				value: Array.from({length: 6}, () => new THREE.Vector3(0, 0, 0)) 
			}, // centreOfArray array with 6 independent Vector3 elements, all (0, 0, 0)
			radius: { value: Array(6).fill(10.0) }, // radius array with 6 elements, all 5.0

			// visible1: { value: true },
			// period1: { value: 0.4 },
			// alpha1: { value: 0.0 },
			// lensletsF1: { value: 1.0 },	// focal length of the individual lenslets
			// additionalF1: { value: 1e10 },	// additional focal length of lenslet array (an additional lens in the same plane)
			// centreOfArray1: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			// radius1: { value: 5.0 },	// radius of array 1
			// visible2: { value: true },
			// period2: { value: 0.4 },
			// alpha2: { value: 0 },
			// lensletsF2: { value: -0.5 },
			// additionalF2: { value: 1e10 },	// additional focal length of lenslet array (an additional lens in the same plane)
			// centreOfArray2: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			// radius2: { value: 5.0 },	// radius of array 2
			idealLenses: { value: true },
			videoFeedUTexture: { value: videoFeedUTexture }, 
			videoFeedETexture: { value: videoFeedETexture }, 
			// cameraLensDistance: { value: cameraLensDistance },
			tanHalfFovHU: { value: 1.0 },
			tanHalfFovVU: { value: 1.0 },
			tanHalfFovHE: { value: 1.0 },
			tanHalfFovVE: { value: 1.0 },
			halfWidthU: { value: 1.0 },
			halfHeightU: { value: 1.0 },
			halfWidthE: { value: 1.0 },
			halfHeightE: { value: 1.0 },
			videoDistance: { value: 1000.0 },	// distance of the image of the video feed from the origin
			focusDistance: { value: 10.0 },
			apertureXHat: { value: new THREE.Vector3(1, 0, 0) },
			apertureYHat: { value: new THREE.Vector3(0, 1, 0) },
			apertureRadius: { value: apertureRadius },
			randomNumbersX: { value: randomNumbersX },
			randomNumbersY: { value: randomNumbersY },
			noOfRays: { value: 1 }
		},
		vertexShader: `
			varying vec3 intersectionPoint;

			void main()	{
				// projectionMatrix, modelViewMatrix, position -> passed in from Three.js
				intersectionPoint = position.xyz;
  				gl_Position = projectionMatrix
					* modelViewMatrix
					* vec4(position, 1.0);
			}
		`,
		fragmentShader: `
			precision highp float;

			varying vec3 intersectionPoint;
			
			// // lenslet array 1
			// uniform bool visible1;
			// uniform float alpha1;	// rotation angle of array 1
			// uniform float period1;	// period of array 1
			// uniform float lensletsF1;	// focal length of array 1
			// uniform float additionalF1;	// additional focal length of lenslet array (an additional lens in the same plane)
			// uniform vec3 centreOfArray1;	// centre of array 1,  and principal point of lenslet (0, 0)
			// uniform float radius1;	// radius of array 1

			// // lenslet array 2
			// uniform bool visible2;
			// uniform float alpha2;	// rotation angle of array 2
			// uniform float period2;	// period of array 2
			// uniform float lensletsF2;	// focal length of array 2
			// uniform float additionalF2;	// additional focal length of lenslet array (an additional lens in the same plane)
			// uniform vec3 centreOfArray2;	// centre of array 2,  and principal point of lenslet (0, 0)
			// uniform float radius2;	// radius of array 2


			// lenslet arrays
			uniform bool visible[6];  // Array to hold visibility flags for both arrays
			uniform float alpha[6];    // Array for rotation angles of both arrays
			uniform float period[6];   // Array for periods of both arrays
			uniform float lensletsF[6];// Array for focal lengths of both arrays
			uniform float additionalF[6]; // Array for additional focal lengths of both arrays
			uniform vec3 centreOfArray[6]; // Array for centers of both arrays
			uniform float radius[6];   // Array for radii of both arrays
			uniform int displayMode;



			uniform bool idealLenses;	// true => use ideal thin lenses; false => use lens holograms

			// video feed from user-facing camera
			uniform sampler2D videoFeedUTexture;
			uniform float halfWidthU;
			uniform float halfHeightU;

			// video feed from environment-facing camera
			uniform sampler2D videoFeedETexture;
			uniform float halfWidthE;
			uniform float halfHeightE;

			// the camera's wide aperture
			uniform float videoDistance;
			uniform float focusDistance;
			uniform int noOfRays;
			uniform vec3 apertureXHat;
			uniform vec3 apertureYHat;
			uniform float apertureRadius;
			uniform float randomNumbersX[100];
			uniform float randomNumbersY[100];
			// uniform float apertureRadius;

			// rotate the 2D vector v by the angle alpha (in radians)
			// from https://gist.github.com/yiwenl/3f804e80d0930e34a0b33359259b556c
			vec2 rotate(vec2 v, float alpha) {
				float s = sin(alpha);
				float c = cos(alpha);
				mat2 m = mat2(c, s, -s, c);
				return m * v;
			}

			// propagate the ray starting at position p and with direction d to the plane z = z0, providing that plane
			// is in the ray's "forward" direction;
			// p becomes the point where the ray intersects p;
			// isForward is set to true or false depending on whether the intersection point is forwards or backwards along the ray
			void propagateForwardToZPlane(
				inout vec3 p, 
				vec3 d, 
				float z0,
				inout bool isForward
			) {
				// calculate the z distance from the ray start position to the array
				float deltaZ = z0 - p.z;

				// is the intersection with the plane in the ray's "forward" direction?
				isForward = (d.z*deltaZ > 0.0);

				// if the intersection is in the forward direction, advance the ray to the plane
				if(isForward) p += d/d.z*deltaZ;	// set p to the intersection point with the plane
			}

			// Calculate the light-ray direction after transmission through a lens or lens hologram.
			// d is the incident light-ray direction;
			// pixy is a 2D vector containing the transverse (x, y) components of the vector I-P,
			// i.e. the vector from the principal point P to the intersection point I;
			// f is the focal length;
			// returns the outgoing light-ray direction
			vec3 lensDeflect(vec3 d, vec2 pixy, float f) {
				if(idealLenses) {
					// ideal thin lens
					// "normalise" the direction such that the magnitude of the z component is 1
					vec3 d1 = d/abs(d.z);

					// the 3D deflected direction comprises the transverse components and a z component of magnitude 1
					// and the same sign as d.z
					return vec3(d1.xy - pixy/f, d1.z);
				} else {
					// lens hologram
					// normalise d
					vec3 dN = d/length(d);
					// transverse components of the outgoing light-ray direction
					vec2 dxy = dN.xy - pixy/f;
	
					// from the transverse direction, construct a 3D vector by setting the z component such that the length
					// of the vector is 1
					return vec3(dxy, sign(d.z)*sqrt(1.0 - dot(dxy, dxy)));
				}
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lens.
			// The (ideal thin) lens, of focal length f, is in a z plane through centreOfLens.
			// It is circular, with the given radius, centred on centreOfLenss.
			void passThroughLens(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfLens, 
				float radius,
				float focalLength
			) {
				bool isForward;
				propagateForwardToZPlane(p, d, centreOfLens.z, isForward);

				if(isForward) {
					// there is an intersection with the plane of this lens in the ray's forward direction

					// does the intersection point lie within the radius?
					vec2 pixy = p.xy - centreOfLens.xy;
					float r2 = dot(pixy, pixy);
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the lens does something to the ray

						// deflect the light-ray direction accordingly and make sure that the sign of the z component remains the same
						lensDeflect(d, pixy, focalLength);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				}
			}

			float findLensletNumber(float u, float uPeriod) {
				return floor(u/uPeriod + 0.5);
			}

			float findLensletCentreCoordinate(float u, float uPeriod, out float lensletNumber) {
				lensletNumber = findLensletNumber(u, uPeriod);
				return uPeriod * lensletNumber;
			}

			// float findLensletCentreCoordinate(float u, float uPeriod) {
			// 	return uPeriod*floor(u/uPeriod+0.5);
			// }

			vec2 rotate2(vec2 v, float cosAlpha, float sinAlpha) {
				float s = (cosAlpha);
				float c = (sinAlpha);
				//float s = sin(cosAlpha);
				//float c = cos(sinAlpha);
				mat2 m = mat2(c, s, -s, c);
				return m * v;
			}

			vec2 findNearestPrincipalPoint(vec2 r, vec2 pp00, float cosAlpha, float sinAlpha, 
				float period1, float period2, out float lensletNumber1, out float lensletNumber2) {
				vec2 rr = rotate2(r, cosAlpha, -sinAlpha);
				vec2 lensletCentreST = vec2(
					findLensletCentreCoordinate(rr.x, period1, lensletNumber1),
					findLensletCentreCoordinate(rr.y, period2, lensletNumber2)
				);
				return rotate2(lensletCentreST, cosAlpha, sinAlpha) + pp00;
			}
			
			vec3 lensletArrayDeflect(vec3 d, vec3 intersectionPoint, inout vec4 b, vec3 principalPoint00, float cosAlpha,
				float sinAlpha, float period, float focalLength, float n) {
				vec2 p00ixy = intersectionPoint.xy - principalPoint00.xy;
				float lensletNumber1;
				float lensletNumber2;
				vec2 pxy = -findNearestPrincipalPoint(p00ixy, principalPoint00.xy, cosAlpha, sinAlpha, period, period, lensletNumber1, lensletNumber2);
				
				if(mod(lensletNumber2, 3.) == n) {
					// lower the brightness factor, giving the light a blue tinge
					b *= vec4(0.9, 0.9, 0.99, 1);
					return lensDeflect(d, intersectionPoint.xy - pxy, focalLength);
				}				
				return d;
			}


			// Find the principal point of the nearest lenslet in a rectangular lenslet array.
			// r is a 2D vector containing the transverse components of the vector from pp00 to the
			// intersection point;
			// pp00 is the principal point of the lenslet with indices (0, 0) -- sort of the array centre;
			// alpha is the angle (in radians) by which the array is rotated
			// period1 and period2 are the periods in the (rotated) x and y directions
			// vec2 findNearestPrincipalPoint(vec2 r, vec2 pp00, float alpha, float period1, float period2) {
			// 	vec2 rr = rotate(r, -alpha);
			// 	vec2 lensletCentreST = vec2(
			// 		findLensletCentreCoordinate(rr.x, period1),
			// 		findLensletCentreCoordinate(rr.y, period2)
			// 	);
			// 	return rotate(lensletCentreST, alpha) + pp00;
			// }

			// vec3 lensletArrayDeflect(vec3 d, vec3 intersectionPoint, vec3 principalPoint00, float alpha, float period, float focalLength) {
			// 	vec2 p00ixy = intersectionPoint.xy - principalPoint00.xy;
			// 	vec2 pxy = findNearestPrincipalPoint(p00ixy, principalPoint00.xy, alpha, period, period);
			// 	// light-ray direction after array
			// 	return lensDeflect(d, intersectionPoint.xy - pxy, focalLength);
			// }

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lenslet array.
			// The lenslet array is in a z plane through centreOfArray; it is simulated as ideal thin lenses.
			// The lenslets, all of focal length f, are arranged in a square array of the given period, 
			// rotated by alpha around the z axis.  The component is circular, with the given radius, centred on centreOfArray.
			void passThroughLensletArray(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfArray, 
				float radius,
				float cosAlpha,
				float sinAlpha,
				float period,
				float f,
				float overallF,
				float n
			) {
				bool isForward;
				propagateForwardToZPlane(p, d, centreOfArray.z, isForward);

				if(isForward) {
					// there is an intersection with the plane of this array in the ray's forward direction
					
					// does the intersection point lie with in the radius?
					vec2 pixy = p.xy - centreOfArray.xy;
					float r2 = dot(pixy, pixy);	// length squared of vector r
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the component does something to the ray

						// deflect the light-ray direction accordingly 

						d = lensletArrayDeflect(d, p, b, centreOfArray, cosAlpha, sinAlpha, period, f, n);
						// d = lensDeflect(d, pixy, overallF);
					} 
				} else b *= vec4(0.99, 0.9, 0.9, 1);	// this shouldn't happen -- give the light a red tinge
			}

			// propagate the ray to the plane of the video feed, which is a z-distance <videoDistance> away,
			// and return either the color of the corresponding video-feed texel or the background color
			vec4 getColorOfVideoFeed(
				inout vec3 p, 
				vec3 d, 
				vec4 b,
				float videoFeedZ,
				sampler2D videoFeedTexture,
				float halfWidth,
				float halfHeight,
				vec4 backgroundColor
			) {
				bool isForward;
				propagateForwardToZPlane(p, d, videoFeedZ, isForward);

				// is the intersection in the ray's forward direction?
				if(isForward) {
					// does the ray intersect the image?
					if((abs(p.x) < halfWidth) && (abs(p.y) < halfHeight))
						// yes, the ray intersects the image; take the pixel colour from the camera's video feed
						return texture2D(videoFeedTexture, vec2(0.5+0.5*p.x/halfWidth, 0.5+0.5*p.y/halfHeight));
					else 
						// the ray doesn't intersect the image
						return backgroundColor;
				}
			}

			void main() {
				// first calculate the point this pixel is focussed on
				vec3 dF = intersectionPoint - cameraPosition;
				vec3 focusPosition = cameraPosition + focusDistance/abs(dF.z)*dF;

				// trace <noOfRays> rays
				gl_FragColor = vec4(0, 0, 0, 0);
				vec4 color;
				for(int i=0; i<noOfRays; i++) {
					// the current ray start position, a random point on the camera's circular aperture
					vec3 p = cameraPosition + apertureRadius*randomNumbersX[i]*apertureXHat + apertureRadius*randomNumbersY[i]*apertureYHat;
	
					// first calculate the current light-ray direction:
					// the ray first passes through focusPosition and then p,
					// so the "backwards" ray direction from the camera to the intersection point is
					//   d = focusPosition - p
					vec3 d = focusPosition - p;
					d = dF.z/d.z*d;
	
					// current brightness factor; this will multiply the colour at the end
					vec4 b = vec4(1.0, 1.0, 1.0, 1.0);
	
					float shift[6];
					if(displayMode == 0) {
						//A模式 第二组和第一组一致
						shift = float[6](0., 0., 1., 1., 2., 2.);
					} else if(displayMode == 1) {
						//B模式 第二组所有整列透镜向后移动一位。
						shift = float[6](0., 2., 1., 0., 2., 1.);
					} else if(displayMode == 2) {
						//C模式 第二组和第一组一致
						shift = float[6](0., 1., 1., 2., 2., 0.);
					}

					if(d.z < 0.0) {
						// the ray is travelling "forwards", in the (-z) direction;
						// pass first through array 1, then array 2, then to environment-facing video feed
						//if(visible[0]) passThroughLensletArray(p, d, b, centreOfArray[0], radius[0], alpha[0], period[0], lensletsF[0], additionalF[0]);
						//if(visible[1]) passThroughLensletArray(p, d, b, centreOfArray[1], radius[1], alpha[1], period[1], lensletsF[1], additionalF[1]);
						
						for (int j = 0; j < 6; j++) {
							if (visible[j]) {
								passThroughLensletArray(p, d, b, centreOfArray[j], radius[j], 
									cos(alpha[j]), sin(alpha[j]), period[j], lensletsF[j], additionalF[j], 
									shift[j]
								);
							}
						}
						
						color = getColorOfVideoFeed(p, d, b, -videoDistance, videoFeedETexture, halfWidthE, halfHeightE, vec4(1, 1, 1, 1.0));	// white
						// if(visible1) passThroughLens(p, d, b, centreOfArray[0],  radius[0], lensletsF[0]);
						// if(visible[1]) passThroughLens(p, d, b, centreOfArray[1],  radius[1], lensletsF[1]);
					} else {
						// the ray is travelling "backwards", in the (+z) direction;
						// pass first through array 2, then array 1, then to user-facing video feed
						//if(visible[1]) passThroughLensletArray(p, d, b, centreOfArray[1], radius[1], alpha[1], period[1], lensletsF[1], additionalF[1]);
						//if(visible[0]) passThroughLensletArray(p, d, b, centreOfArray[0], radius[0], alpha[0], period[0], lensletsF[0], additionalF[0]);
						
						for (int j = 5; j >= 0; j--) {
							if (visible[j]) {
								passThroughLensletArray(p, d, b, centreOfArray[j], radius[j], cos(alpha[j]), 
									sin(alpha[j]), period[j], lensletsF[j], additionalF[j],
									shift[j]
								);
							}
						}
						
						color = getColorOfVideoFeed(p, d, b, videoDistance, videoFeedUTexture, halfWidthU, halfHeightU, vec4(1, 0, 0, 1.0));	// white
					}
		
					// finally, multiply by the brightness factor and add to gl_FragColor
					gl_FragColor += b*color;
				}
					
				gl_FragColor /= float(noOfRays);
			}
		`
	});
	raytracingSphere = new THREE.Mesh( geometry, raytracingSphereShaderMaterial ); 
	scene.add( raytracingSphere );
}

function addEventListenersEtc() {
	// handle device orientation
	// window.addEventListener("deviceorientation", handleOrientation, true);

	// handle window resize
	window.addEventListener("resize", onWindowResize, false);

	// handle screen-orientation (landscape/portrait) change
	screen.orientation.addEventListener( "change", recreateVideoFeeds );

	// share button functionality
	document.getElementById('takePhotoButton').addEventListener('click', takePhoto);

	// toggle fullscreen button functionality
	document.getElementById('fullscreenButton').addEventListener('click', toggleFullscreen);

	// info button functionality
	document.getElementById('infoButton').addEventListener('click', toggleInfoVisibility);

	// back button functionality
	document.getElementById('backButton').addEventListener('click', showLivePhoto);
	document.getElementById('backButton').style.visibility = "hidden";

	// share button
	document.getElementById('shareButton').addEventListener('click', share);
	document.getElementById('shareButton').style.visibility = "hidden";
	if(!(navigator.share)) document.getElementById('shareButton').src="./shareButtonUnavailable.png";
	// if(!(navigator.share)) document.getElementById('shareButton').style.opacity = 0.3;

	// delete button
	document.getElementById('deleteButton').addEventListener('click', deleteStoredPhoto);
	document.getElementById('deleteButton').style.visibility = "hidden";

	// hide the thumbnail for the moment
	document.getElementById('storedPhotoThumbnail').addEventListener('click', showStoredPhoto);
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('storedPhoto').addEventListener('click', showLivePhoto);
	document.getElementById('storedPhoto').style.visibility = "hidden";
	// showingStoredPhoto = false;
}

// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_animation_skinning_additive_blending.html
function createGUI() {
	// const 
	gui = new GUI();
	const paramsMode = {
		mode: 'A'  // 默认模式A
	};
	gui.add(paramsMode, 'mode', ['A', 'B', 'C']).onChange(value => {
		if (value === 'A') {
			raytracingSphereShaderMaterial.uniforms.displayMode.value = 0;
		} else if (value === 'B') {
			raytracingSphereShaderMaterial.uniforms.displayMode.value = 1;
		} else if (value === 'C') {
			raytracingSphereShaderMaterial.uniforms.displayMode.value = 2;
		}
		updateUniforms();
		
	});
	
	// gui.hide();

	// const params1 = {
	// 	'Visible': raytracingSphereShaderMaterial.uniforms.visible1.value,
	// 	'Lenslets <i>f</i><sub>1</sub>': raytracingSphereShaderMaterial.uniforms.lensletsF1.value,
	// 	'tan<sup>-1</sup>(additional <i>F</i><sub>1</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF1.value),
	// 	'Period, <i>p</i><sub>1</sub>': raytracingSphereShaderMaterial.uniforms.period1.value,
	// 	'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha1.value / Math.PI * 180.
	// };
	// const params2 = {
	// 	'Visible': raytracingSphereShaderMaterial.uniforms.visible2.value,
	// 	'Lenslets <i>f</i><sub>2</sub>': raytracingSphereShaderMaterial.uniforms.lensletsF2.value,
	// 	'tan<sup>-1</sup>(additional <i>F</i><sub>2</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF2.value),
	// 	'&Delta;<i>p</i> (<i>p</i><sub>2</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)': deltaPeriod,
	// 	'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha2.value / Math.PI * 180.,
	// 	'Offset from confocal': offsetFromConfocal
	// };
	// const params3 = {
	// 	'Visible': raytracingSphereShaderMaterial.uniforms.visible2.value,
	// 	'Lenslets <i>f</i><sub>3</sub>': raytracingSphereShaderMaterial.uniforms.lensletsF2.value,
	// 	'tan<sup>-1</sup>(additional <i>F</i><sub>3</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF2.value),
	// 	'&Delta;<i>p</i> (<i>p</i><sub>3</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)': deltaPeriod,
	// 	'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha2.value / Math.PI * 180.,
	// 	'Offset from confocal': offsetFromConfocal
	// };
	// const params4 = {
	// 	'Visible': raytracingSphereShaderMaterial.uniforms.visible2.value,
	// 	'Lenslets <i>f</i><sub>4</sub>': raytracingSphereShaderMaterial.uniforms.lensletsF2.value,
	// 	'tan<sup>-1</sup>(additional <i>F</i><sub>4</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF2.value),
	// 	'&Delta;<i>p</i> (<i>p</i><sub>4</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)': deltaPeriod,
	// 	'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha2.value / Math.PI * 180.,
	// 	'Offset from confocal': offsetFromConfocal
	// };
	// const params5 = {
	// 	'Visible': raytracingSphereShaderMaterial.uniforms.visible2.value,
	// 	'Lenslets <i>f</i><sub>5</sub>': raytracingSphereShaderMaterial.uniforms.lensletsF2.value,
	// 	'tan<sup>-1</sup>(additional <i>F</i><sub>5</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF2.value),
	// 	'&Delta;<i>p</i> (<i>p</i><sub>5</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)': deltaPeriod,
	// 	'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha2.value / Math.PI * 180.,
	// 	'Offset from confocal': offsetFromConfocal
	// };
	// const params6 = {
	// 	'Visible': raytracingSphereShaderMaterial.uniforms.visible2.value,
	// 	'Lenslets <i>f</i><sub>6</sub>': raytracingSphereShaderMaterial.uniforms.lensletsF2.value,
	// 	'tan<sup>-1</sup>(additional <i>F</i><sub>6</sub>)': Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF2.value),
	// 	'&Delta;<i>p</i> (<i>p</i><sub>6</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)': deltaPeriod,
	// 	'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha2.value / Math.PI * 180.,
	// 	'Offset from confocal': offsetFromConfocal
	// };

	const params = {
		// 'Swap arrays': swapArrays,
		'Lenslet type': 'Ideal thin',
		'Ideal lenses': raytracingSphereShaderMaterial.uniforms.idealLenses.value,
		'Horiz. FOV (&deg;)': fovScreen,
		'Aperture radius': apertureRadius,
		'tan<sup>-1</sup>(focus. dist.)': Math.atan(focusDistance),
		'No of rays': noOfRays,
		'Env.-facing cam. (&deg;)': fovVideoFeedE,
		'User-facing cam. (&deg;)': fovVideoFeedU,
		'tan<sup>-1</sup>(video dist.)': Math.atan(raytracingSphereShaderMaterial.uniforms.videoDistance.value),
		'Point (virtual) cam. forward (in -<b>z</b> direction)': pointForward,
		'Show/hide info': toggleInfoVisibility,
		'Restart video streams': function() { 
			recreateVideoFeeds(); 
			postStatus("Restarting video stream");
		}
	}

	// const folderArray1 = gui.addFolder( 'Lenslet array 1 (near; <i>z</i><sub>1</sub> = 0)' );

	// folderArray1.add( params1, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible1.value = v; } );
	// folderArray1.add( params1, 'Lenslets <i>f</i><sub>1</sub>', -1, 1).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.lensletsF1.value = f; } );
	// folderArray1.add( params1, 'tan<sup>-1</sup>(additional <i>F</i><sub>1</sub>)', -0.5*Math.PI, 0.5*Math.PI).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.additionalF1.value = Math.tan(f); } );
	// folderArray1.add( params1, 'Period, <i>p</i><sub>1</sub>', 0.01, 1).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.period1.value = p; } );
	// folderArray1.add( params1, 'Rotation angle (&deg;)', -10, 10).onChange( (alpha) => { raytracingSphereShaderMaterial.uniforms.alpha1.value = alpha/180.0*Math.PI; } );

	// const folderArray2 = gui.addFolder( 'Lenslet array 2 (far)' );

	// folderArray2.add( params2, 'Visible').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible2.value = v; } );
	// folderArray2.add( params2, 'Lenslets <i>f</i><sub>2</sub>', -1, 1).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.lensletsF2.value = f; } );
	// folderArray2.add( params2, 'tan<sup>-1</sup>(additional <i>F</i><sub>2</sub>)', -0.5*Math.PI, 0.5*Math.PI).onChange( (f) => { raytracingSphereShaderMaterial.uniforms.additionalF2.value = Math.tan(f); } );
	// folderArray2.add( params2, '&Delta;<i>p</i> (<i>p</i><sub>2</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)', -0.1, 0.1).onChange( (p) => { deltaPeriod = p; } );
	// folderArray2.add( params2, 'Rotation angle (&deg;)', -10, 10).onChange( (alpha) => { raytracingSphereShaderMaterial.uniforms.alpha2.value = alpha/180.0*Math.PI; } );
	// folderArray2.add( params2, 'Offset from confocal', -0.25, 0.25).onChange( (o) => { offsetFromConfocal = o; } );


	const numArrays = 6;
	//const deltaPeriod = 0.1; // 假设这个变量已经定义
	//const offsetFromConfocal = 0.2; // 假设这个变量已经定义

	// 创建数组
	const paramsArray = [];
	for (let i = 0; i < numArrays; i++) {
		paramsArray[i] = {
			'Visible': raytracingSphereShaderMaterial.uniforms.visible.value[i],
			[`Lenslets <i>f</i><sub>${i + 1}</sub>`]: raytracingSphereShaderMaterial.uniforms.lensletsF.value[i],
			[`tan<sup>-1</sup>(additional <i>F</i><sub>${i + 1}</sub>)`]: Math.atan(raytracingSphereShaderMaterial.uniforms.additionalF.value[i]),
			[`&Delta;<i>p</i> (<i>p</i><sub>${i + 1}</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)`]: deltaPeriod,
			'Rotation angle (&deg;)': raytracingSphereShaderMaterial.uniforms.alpha.value[i] / Math.PI * 180.,
			'Offset from confocal': offsetFromConfocal
		};
	}

	// 创建文件夹和控件
	for (let i = 0; i < numArrays; i++) {
		const folderArray = gui.addFolder(`Lenslet array ${i + 1} (far)`);

		folderArray.add(paramsArray[i], 'Visible').onChange((v) => { raytracingSphereShaderMaterial.uniforms.visible.value[i] = v; });
		folderArray.add(paramsArray[i], `Lenslets <i>f</i><sub>${i + 1}</sub>`, -10, 10).onChange((f) => { raytracingSphereShaderMaterial.uniforms.lensletsF.value[i] = f; });
		// folderArray.add(paramsArray[i], `tan<sup>-1</sup>(additional <i>F</i><sub>${i + 1}</sub>)`, -0.5 * Math.PI, 0.5 * Math.PI).onChange((f) => { raytracingSphereShaderMaterial.uniforms.additionalF.value[i] = Math.tan(f); });
		// folderArray.add(paramsArray[i], `&Delta;<i>p</i> (<i>p</i><sub>${i + 1}</sub> = <i>p</i><sub>1</sub> + &Delta;<i>p</i>)`, -0.1, 0.1).onChange((p) => { deltaPeriod = p; });
		// folderArray.add(paramsArray[i], 'Rotation angle (&deg;)', -10, 10).onChange((alpha) => { raytracingSphereShaderMaterial.uniforms.alpha.value[i] = alpha / 180.0 * Math.PI; });
		// folderArray.add(paramsArray[i], 'Offset from confocal', -0.25, 0.25).onChange((o) => { offsetFromConfocal = o; });
	}

	

	const folderDevice = gui.addFolder( 'Device cameras horiz. FOV' );
	folderDevice.add( params, 'Env.-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedE = fov; updateUniforms(); });   
	folderDevice.add( params, 'User-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedU = fov; updateUniforms(); });   
	folderDevice.close();

	const folderVirtualCamera = gui.addFolder( 'Virtual camera' );
	folderVirtualCamera.add( params, 'Horiz. FOV (&deg;)', 10, 170, 1).onChange( setScreenFOV );
	folderVirtualCamera.add( params, 'Aperture radius', 0.0, 1.0).onChange( (r) => { apertureRadius = r; } );
	folderVirtualCamera.add( params, 'tan<sup>-1</sup>(focus. dist.)', 
		//Math.atan(0.1), 
		-0.5*Math.PI,
		0.5*Math.PI
	).onChange( (a) => { focusDistance = Math.tan(a); } );
	folderVirtualCamera.add( params, 'No of rays', 1, 100, 1).onChange( (n) => { noOfRays = n; } );

	const folderSettings = gui.addFolder( 'Other controls' );
	folderSettings.add( params, 'tan<sup>-1</sup>(video dist.)', Math.atan(0.1), 0.5*Math.PI).onChange( (a) => { raytracingSphereShaderMaterial.uniforms.videoDistance.value = Math.tan(a); } );
	folderSettings.add( params, 'Lenslet type', { 'Ideal thin': true, 'Phase hologram': false } ).onChange( (t) => { raytracingSphereShaderMaterial.uniforms.idealLenses.value = t; });
	// folderSettings.add( params, 'Ideal lenses').onChange( (b) => { raytracingSphereShaderMaterial.uniforms.idealLenses.value = b; } );
	folderSettings.add( params, 'Point (virtual) cam. forward (in -<b>z</b> direction)');
	folderSettings.add( params, 'Show/hide info');
	folderSettings.add( params, 'Restart video streams');
	folderSettings.close();

	const distanceSettings = gui.addFolder( 'Distance controls' );
	const paramsDistance = {
		'Distance between lens 1 & 2': distance12,  // 1号和2号之间的距离
		'Distance between lens 3 & 4': distance34,  // 3号和4号之间的距离
		'Distance between lens 5 & 6': distance56,  // 5号和6号之间的距离
		'Gap between lens 2 & 3': gapBetween2And3,       // 2号和3号之间的间隔
		'Gap between lens 4 & 5': gapBetween4And5,       // 4号和5号之间的间隔
		'Unified rotation angle': 0           // 统一旋转角度
	};
	
	
	
	// 控制1号与2号、3号与4号、5号与6号之间的距离
	distanceSettings.add(paramsDistance, 'Distance between lens 1 & 2', 0, 10).onChange(value => {
		distance12 = value;
		updateUniforms();  // 每次改变时调用 updateUniforms 以更新透镜的参数
	});
	distanceSettings.add(paramsDistance, 'Distance between lens 3 & 4', 0, 10).onChange(value => {
		distance34 = value;
		updateUniforms();
	});
	distanceSettings.add(paramsDistance, 'Distance between lens 5 & 6', 0, 10).onChange(value => {
		distance56 = value;
		updateUniforms();
	});
	
	// 控制2号与3号、4号与5号之间的间隔
	distanceSettings.add(paramsDistance, 'Gap between lens 2 & 3', 0, 5).onChange(value => {
		gapBetween2And3 = value;
		updateUniforms();
	});
	distanceSettings.add(paramsDistance, 'Gap between lens 4 & 5', 0, 5).onChange(value => {
		gapBetween4And5 = value;
		updateUniforms();
	});
	
	// 控制所有透镜的统一旋转角度
	distanceSettings.add(paramsDistance, 'Unified rotation angle', -Math.PI, Math.PI).onChange(value => {
		rotationAngle = value;
		updateUniforms();
	});


}

/**
 * @param {*} fov	The larger of the camera's horizontal and vertical FOV, in degrees
 * 
 * Set the larger FOV of the screen/window to fov.
 * 
 * Depending on the screen/window's FOV, fov is either the horizontal fov (if screen width > screen height)
 * or the vertical fov (if screen width < screen height).
 */
function setScreenFOV(fov) {
	fovScreen = fov;

	screenChanged();
}

// function swapArrays() {
// 	const visible3 = lensletArrayShaderMaterial.uniforms.visible1.value;
// 	const focalLength3 = lensletArrayShaderMaterial.uniforms.lensletsF1.value;
// 	const period3 = lensletArrayShaderMaterial.uniforms.period1.value;
// 	const alpha3 = lensletArrayShaderMaterial.uniforms.alpha1.value;

// 	lensletArrayShaderMaterial.uniforms.visible1.value = lensletArrayShaderMaterial.uniforms.visible2.value;
// 	lensletArrayShaderMaterial.uniforms.lensletsF1.value = lensletArrayShaderMaterial.uniforms.lensletsF2.value;
// 	lensletArrayShaderMaterial.uniforms.period1.value = lensletArrayShaderMaterial.uniforms.period2.value;
// 	lensletArrayShaderMaterial.uniforms.alpha1.value = lensletArrayShaderMaterial.uniforms.alpha2.value;

// 	lensletArrayShaderMaterial.uniforms.visible2.value = visible3;
// 	lensletArrayShaderMaterial.uniforms.lensletsF2.value = focalLength3;
// 	lensletArrayShaderMaterial.uniforms.period2.value = period3;
// 	lensletArrayShaderMaterial.uniforms.alpha2.value = alpha3;
// }

/** 
 * Reset the aspect ratio and FOV of the virtual cameras.
 * 
 * Call if the window size has changed (which also happens when the screen orientation changes)
 * or if camera's FOV has changed
 */
function screenChanged() {
	// alert(`new window size ${window.innerWidth} x ${window.innerHeight}`);

	// in case the screen size has changed
	if(renderer) renderer.setSize(window.innerWidth, window.innerHeight);

	// if the screen orientation changes, width and height swap places, so the aspect ratio changes
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera.aspect = windowAspectRatio;

	// fovS is the screen's horizontal or vertical FOV, whichever is greater;
	// re-calculate the camera FOV, which is the *vertical* fov
	let verticalFOV;
	if(windowAspectRatio > 1.0) {
		// fovS is horizontal FOV; convert to get correct vertical FOV
		verticalFOV = 2.0*Math.atan(Math.tan(0.5*fovScreen*Math.PI/180.0)/windowAspectRatio)*180.0/Math.PI;
	} else {
		// fovS is already vertical FOV
		verticalFOV = fovScreen;
	}
	camera.fov = verticalFOV;

	// make sure the camera changes take effect
	camera.updateProjectionMatrix();
}

function  pointForward() {
	let r = camera.position.length();
	camera.position.x = 0;
	camera.position.y = 0;
	camera.position.z = r;
	controls.update();
	postStatus('Pointing camera forwards (in -<b>z</b> direction)');
}

function onWindowResize() {
	screenChanged();
	postStatus(`window size ${window.innerWidth} &times; ${window.innerHeight}`);	// debug
}

// // see https://developer.mozilla.org/en-US/docs/Web/API/ScreenOrientation/change_event
function recreateVideoFeeds() {
	// stop current video streams...
	videoFeedE.srcObject.getTracks().forEach(function(track) { track.stop(); });
	videoFeedU.srcObject.getTracks().forEach(function(track) { track.stop(); });

	// ... and re-create new ones, hopefully of the appropriate size
	createVideoFeeds();
}

function addOrbitControls() {
	// controls

	controls = new OrbitControls( camera, renderer.domElement );
	// controls = new OrbitControls( cameraOutside, renderer.domElement );
	controls.listenToKeyEvents( window ); // optional

	//controls.addEventListener( 'change', render ); // call this only in static scenes (i.e., if there is no animation loop)
	controls.addEventListener( 'change', cameraPositionChanged );

	controls.enableDamping = false; // an animation loop is required when either damping or auto-rotation are enabled
	controls.dampingFactor = 0.05;

	controls.enablePan = true;
	controls.enableZoom = true;

	controls.maxPolarAngle = Math.PI;
}

function cameraPositionChanged() {
	postStatus(`Camera position (${camera.position.x.toPrecision(2)}, ${camera.position.y.toPrecision(2)}, ${camera.position.z.toPrecision(2)})`);
	// counter = 0;
	// keep the raytracing sphere centred on the camera position
	// raytracingSphere.position.copy(camera.position.clone());	// TODO this doesn't seem to work as intended!?
}

async function toggleFullscreen() {
	if (!document.fullscreenElement) {
		document.documentElement.requestFullscreen().catch((err) => {
			postStatus(
				`Error attempting to enable fullscreen mode: ${err.message} (${err.name})`,
			);
		});
		// allow screen orientation changes
		// screen.orientation.unlock();
	} else {
		document.exitFullscreen();
	}
}

function showStoredPhoto() {
	gui.hide();
	renderer.domElement.style.visibility = "hidden";
	document.getElementById('takePhotoButton').style.visibility = "hidden";
	// document.getElementById('changePositionButton').style.visibility = "hidden";
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('backButton').style.visibility = "visible";
	document.getElementById('shareButton').style.visibility = "visible";
	document.getElementById('deleteButton').style.visibility = "visible";
	document.getElementById('storedPhoto').style.visibility = "visible";
	showingStoredPhoto = true;

	postStatus('Showing stored photo, '+storedPhotoDescription);
}

function showLivePhoto() {
	gui.show();
	renderer.domElement.style.visibility = "visible";
	document.getElementById('takePhotoButton').style.visibility = "visible";
	// document.getElementById('changePositionButton').style.visibility = "visible";
	if(storedPhoto) document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	document.getElementById('backButton').style.visibility = "hidden";
	document.getElementById('shareButton').style.visibility = "hidden";
	document.getElementById('deleteButton').style.visibility = "hidden";
	document.getElementById('storedPhoto').style.visibility = "hidden";
	showingStoredPhoto = false;

	postStatus('Showing live image');
}

function deleteStoredPhoto() {
	storedPhoto = null;

	showLivePhoto();

	postStatus('Stored photo deleted; showing live image');
}

function takePhoto() {
	try {
		click.play();

		storedPhoto = renderer.domElement.toDataURL('image/png');
		storedPhotoInfoString = getInfoString();

		storedPhotoDescription = name;
		// 
		document.getElementById('storedPhoto').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	
		postStatus('Photo taken; click thumbnail to view and share');
	} catch (error) {
		console.error('Error:', error);
	}	
}

async function share() {
	try {
		fetch(storedPhoto)
		.then(response => response.blob())
		.then(blob => {
			const file = new File([blob], name+storedPhotoDescription+'.png', { type: blob.type });

			// Use the Web Share API to share the screenshot
			if (navigator.share) {
				navigator.share({
					title: storedPhotoDescription,
					files: [file],
				});
			} else {
				postStatus('Sharing is not supported by this browser.');
			}	
		})
		.catch(error => {
			console.error('Error:', error);
			postStatus(`Error: ${error}`);
		});
	} catch (error) {
		console.error('Error:', error);
	}
}

/** 
 * Add a text field to the bottom left corner of the screen
 */
function createStatus() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	status.style.position = 'absolute';
	status.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	status.style.color = "White";
	status.style.fontFamily = "Arial";
	status.style.fontSize = "9pt";
	postStatus("Welcome!");
	status.style.bottom = 0 + 'px';
	status.style.left = 0 + 'px';
	status.style.zIndex = 1;
	document.body.appendChild(status);	
}

function postStatus(text) {
	status.innerHTML = '&nbsp;'+text;
	console.log('status: '+text);

	// show the text only for 3 seconds
	statusTime = new Date().getTime();
	setTimeout( () => { if(new Date().getTime() - statusTime > 2999) status.innerHTML = '&nbsp;'+name+', University of Glasgow, <a href="https://github.com/jkcuk/'+name+'">https://github.com/jkcuk/'+name+'</a>' }, 3000);
}

function getInfoString() {
	return "not implemented";
	return `Lenslet array 1 (the closer array, when seen in "forward" direction)<br>` +
		`&nbsp;&nbsp;Visible `+ (raytracingSphereShaderMaterial.uniforms.visible1.value?'&check;':'&cross;')+`<br>` +
		`&nbsp;&nbsp;Period = ${raytracingSphereShaderMaterial.uniforms.period1.value.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Rotation angle = ${raytracingSphereShaderMaterial.uniforms.alpha1.value.toPrecision(4)}&deg;<br>` +
		`&nbsp;&nbsp;Focal length = ${raytracingSphereShaderMaterial.uniforms.lensletsF1.value.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Radius = ${raytracingSphereShaderMaterial.uniforms.radius1.value.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Centre of array = (${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.z.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;Focal length of additional lens in same plane = ${raytracingSphereShaderMaterial.uniforms.additionalF1.value.toPrecision(4)}<br>` +		
		`Lenslet array 2 (the farther array, when seen in "forward" direction)<br>` +
		`&nbsp;&nbsp;Visible `+ (raytracingSphereShaderMaterial.uniforms.visible2.value?'&check;':'&cross;')+`<br>` +
		`&nbsp;&nbsp;Period = ${raytracingSphereShaderMaterial.uniforms.period2.value.toPrecision(4)} (&Delta;<i>p</i> = ${deltaPeriod.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;Rotation angle = ${raytracingSphereShaderMaterial.uniforms.alpha2.value.toPrecision(4)}&deg;<br>` +
		`&nbsp;&nbsp;Focal length = ${raytracingSphereShaderMaterial.uniforms.lensletsF2.value.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Radius = ${raytracingSphereShaderMaterial.uniforms.radius2.value.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Centre of array = (${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.z.toPrecision(4)}) (offset from confocal = ${offsetFromConfocal.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;Focal length of additional lens in same plane = ${raytracingSphereShaderMaterial.uniforms.additionalF2.value.toPrecision(4)}<br>` +		
		'Lenslet type: '+(raytracingSphereShaderMaterial.uniforms.idealLenses.value?'Ideal thin lenses':'Phase holograms') + "<br>" +
		`Video feeds<br>` +
		`&nbsp;&nbsp;Distance from origin = ${raytracingSphereShaderMaterial.uniforms.videoDistance.value.toPrecision(4)}<br>` +	// (user-facing) camera
		`&nbsp;&nbsp;Horizontal fields of view (when seen from the origin)<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;User-facing camera = ${fovVideoFeedU.toPrecision(4)}&deg;<br>` +	// (user-facing) camera
		`&nbsp;&nbsp;&nbsp;&nbsp;Environment-facing camera = ${fovVideoFeedE.toPrecision(4)}&deg;<br>` +	// (environment-facing) camera
		`Virtual camera<br>` +
		`&nbsp;&nbsp;Position = (${camera.position.x.toPrecision(4)}, ${camera.position.y.toPrecision(4)}, ${camera.position.z.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;Horiz. FOV = ${fovScreen.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Aperture radius = ${apertureRadius.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Focussing distance = ${focusDistance.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Number of rays = ${noOfRays}`
		// `apertureXHat = (${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.z.toPrecision(4)})<br>` +
		// `apertureYHat = (${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.z.toPrecision(4)})`
		;
		console.log("*");
}

function refreshInfo() {
	if(showingStoredPhoto) setInfo( storedPhotoInfoString );
	else setInfo( getInfoString() );

	if(info.style.visibility == "visible") setTimeout( refreshInfo , 100);	// refresh again a while
}

/** 
 * Add a text field to the top left corner of the screen
 */
function createInfo() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	info.style.position = 'absolute';
	info.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	info.style.color = "White";
	info.style.fontFamily = "Arial";
	info.style.fontSize = "9pt";
	info.innerHTML = "-- nothing to show (yet) --";
	info.style.top = 60 + 'px';
	info.style.left = 0 + 'px';
	info.style.zIndex = 1;
	document.body.appendChild(info);
	info.style.visibility = "hidden";
}

function setInfo(text) {
	info.innerHTML = text;
	console.log('info: '+text);
	// // show the text only for 3 seconds
	// infoTime = new Date().getTime();
	// setTimeout( () => { if(new Date().getTime() - infoTime > 2999) info.innerHTML = `` }, 3000);
	// info.style.visibility = "visible";
}

function toggleInfoVisibility() {
	switch(info.style.visibility) {
		case "visible":
			info.style.visibility = "hidden";
			break;
		case "hidden":
		default:
			info.style.visibility = "visible";
			refreshInfo();
	}
}