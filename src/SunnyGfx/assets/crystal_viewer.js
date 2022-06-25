'use strict';

(async function () {
    const THREE = await import("https://cdn.skypack.dev/three@0.132.2");
    const OrbitControls = (await import("https://cdn.skypack.dev/three@0.132.2/examples/jsm/controls/OrbitControls.js")).OrbitControls;

    // In development, a default key is useful.
    let key = "0";
    // In production, this line will be replaced to define `key` as a unique
    // string identifier. This uniqueness is necessary to support multiple
    // visualizations within the same web page.
    'DEFINE_KEY';

    // Get a reference to the container element that will hold our scene
    const container = document.querySelector(`#scene-container-${key}`);

    function showError(err) {
        container.innerHTML = `<div class="error" style="color:red; font-size:20px;">${err}</div>`;
        throw err;
    }

    try {
        // create a Scene
        const scene = new THREE.Scene();
        scene.background = new THREE.Color("#EEEEEE");

        // read in data from Sunny 
        const data = JSON.parse(document.getElementById(`data-${key}`).textContent);
        const cellTypes = data.cellTypes;
        const bondColors = data.bondColors;
        const bondLabels = data.bondLabels;
        const bondTypeIds = data.bondTypeIds;
        const bondVecs = data.bondVecs;
        const lattVecs = data.lattVecs;
        const basisVecs = data.basisVecs;
        const lattCells = data.lattCells;
        const atomsPerCell = data.atomsPerCell;

        // CPKM colors for all atoms
        var atomColors = new Object();
        atomColors["H"] = "#FFFFFF";  atomColors["He"] = "#D9FFFF"; atomColors["Li"] = "#CC80FF"; atomColors["Be"] = "#C2FF00";
        atomColors["B"] = "#FFB5B5";  atomColors["C"] = "#909090";  atomColors["N"] = "#3050F8";  atomColors["O"] = "#FF0D0D";
        atomColors["F"] = "#90E050";  atomColors["Ne"] = "#B3E3F5"; atomColors["Na"] = "#AB5CF2"; atomColors["Mg"] = "#8AFF00";
        atomColors["Al"] = "#BFA6A6"; atomColors["Si"] = "#F0C8A0"; atomColors["P"] = "#FF8000";  atomColors["S"] = "#FFFF30";
        atomColors["Cl"] = "#1FF01F"; atomColors["Ar"] = "#80D1E3"; atomColors["K"] = "#8F40D4";  atomColors["Ca"] = "#3DFF00";
        atomColors["Sc"] = "#E6E6E6"; atomColors["Ti"] = "#BFC2C7"; atomColors["V"] = "#A6A6AB";  atomColors["Cr"] = "#8A99C7";
        atomColors["Mn"] = "#9C7AC7"; atomColors["Fe"] = "#E06633"; atomColors["Co"] = "#F090A0"; atomColors["Ni"] = "#50D050";
        atomColors["Cu"] = "#C88033"; atomColors["Zn"] = "#7D80B0"; atomColors["Ga"] = "#C28F8F"; atomColors["Ge"] = "#668F8F";
        atomColors["As"] = "#BD80E3"; atomColors["Se"] = "#FFA100"; atomColors["Br"] = "#A62929"; atomColors["Kr"] = "#5CB8D1";
        atomColors["Rb"] = "#702EB0"; atomColors["Sr"] = "#00FF00"; atomColors["Y"] = "#94FFFF";  atomColors["Zr"] = "#94E0E0";
        atomColors["Nb"] = "#73C2C9"; atomColors["Mo"] = "#54B5B5"; atomColors["Tc"] = "#3B9E9E"; atomColors["Ru"] = "#248F8F";
        atomColors["Rh"] = "#0A7D8C"; atomColors["Pd"] = "#006985"; atomColors["Ag"] = "#C0C0C0"; atomColors["Cd"] = "#FFD98F";
        atomColors["In"] = "#A67573"; atomColors["Sn"] = "#668080"; atomColors["Sb"] = "#9E63B5"; atomColors["Te"] = "#D47A00";
        atomColors["I"] = "#940094";  atomColors["Xe"] = "#429EB0"; atomColors["Cs"] = "#57178F"; atomColors["Ba"] = "#00C900";
        atomColors["La"] = "#70D4FF"; atomColors["Ce"] = "#FFFFC7"; atomColors["Pr"] = "#D9FFC7"; atomColors["Nd"] = "#C7FFC7";
        atomColors["Pm"] = "#A3FFC7"; atomColors["Sm"] = "#8FFFC7"; atomColors["Eu"] = "#61FFC7"; atomColors["Gd"] = "#45FFC7";
        atomColors["Tb"] = "#30FFC7"; atomColors["Dy"] = "#1FFFC7"; atomColors["Ho"] = "#00FF9C"; atomColors["Er"] = "#00E675";
        atomColors["Tm"] = "#00D452"; atomColors["Yb"] = "#00BF38"; atomColors["Lu"] = "#00AB24"; atomColors["Hf"] = "#4DC2FF";
        atomColors["Ta"] = "#4DA6FF"; atomColors["W"] = "#2194D6";  atomColors["Re"] = "#267DAB"; atomColors["Os"] = "#266696";
        atomColors["Ir"] = "#175487"; atomColors["Pt"] = "#D0D0E0"; atomColors["Au"] = "#FFD123"; atomColors["Hg"] = "#B8B8D0";
        atomColors["Tl"] = "#A6544D"; atomColors["Pb"] = "#575961"; atomColors["Bi"] = "#9E4FB5"; atomColors["Po"] = "#AB5C00";
        atomColors["At"] = "#754F45"; atomColors["Rn"] = "#428296"; atomColors["Fr"] = "#420066"; atomColors["Ra"] = "#007D00";
        atomColors["Ac"] = "#70ABFA"; atomColors["Th"] = "#00BAFF"; atomColors["Pa"] = "#00A1FF"; atomColors["U"] = "#008FFF";
        atomColors["Np"] = "#0080FF"; atomColors["Pu"] = "#006BFF"; atomColors["Am"] = "#545CF2"; atomColors["Cm"] = "#785CE3";
        atomColors["Bk"] = "#8A4FE3"; atomColors["Cf"] = "#A136D4"; atomColors["Es"] = "#B31FD4"; atomColors["Fm"] = "#B31FBA";
        atomColors["Md"] = "#B30DA6"; atomColors["No"] = "#BD0D87"; atomColors["Lr"] = "#C70066"; atomColors["Rf"] = "#CC0059";
        atomColors["Db"] = "#D1004F"; atomColors["Sg"] = "#D90045"; atomColors["Bh"] = "#E00038"; atomColors["Hs"] = "#E6002E";
        atomColors["Mt"] = "#EB0026";

        // get magnitude of vector
        function mag(v){
            var s = 0.0;
            for(let i=0; i < v.length; i++){
                s += v[i]**2;
            }
            return Math.sqrt(s);
        }
        
        // vector difference
        function diff(vi, vf){
            var d = [0, 0, 0];
            for(let i=0; i < 3; i++){
                d[i] = vf[i] - vi[i];
            }
            return d;
        }

        // get unique list of atoms types
        function onlyUnique(value, index, self) {
            return self.indexOf(value) === index;
        }
        var types = cellTypes.filter(onlyUnique);

        // strings for bond types
        var bondTypes = [];
        for(let i=0; i < bondTypeIds.length; i++){
            bondTypes.push([cellTypes[bondTypeIds[i][0]], cellTypes[bondTypeIds[i][1]]]);
        }

        // coordinate of central unit cell 
        var centerCell = [0,0,0];
        for(let i=0; i < 3; i++){
            for(let j=0; j < 3; j++){
                centerCell[j] += Math.floor((lattCells[i])/2) * lattVecs[i][j];
            }
        }
        // show bonds for sublattices at these positions
        var bondCenters = [];
        for(let a=0; a < cellTypes.length; a++){
            bondCenters.push(centerCell.map((v,d) => v+basisVecs[a][d]));
        }

        // get max lattice spacing
        var maxLattUnit = Math.max(...lattVecs.map(mag));

        // get max bond length
        var maxBondMag = Math.max(...bondVecs.map((v,i) => mag(v[0])));

        // radial function for setting atom transparency
        function alpha(v, c){
            if(maxBondMag < maxLattUnit){
                return 1.0;
            }
            var dr = mag(diff(v, c));
            return (dr > maxBondMag) ? 1.0/(dr-maxBondMag) : 1.0;
        }

        // add atoms to scene as spheres
        var atoms = [];
        var typeCounts = [];
        for(let i=0; i < types.length; i++){
            atoms.push([]);
            typeCounts.push(0);
        }

        var coordMinMax = [ [Infinity, -Infinity], [Infinity, -Infinity], [Infinity, -Infinity] ];
        var cellPos = [0.0, 0.0, 0.0];

        var atomRadius = 0.075 * maxLattUnit;

        for(let i=0; i < lattCells[0]; i++){
            for(let j=0; j < lattCells[1]; j++){
                for(let k=0; k < lattCells[2]; k++){

                    for(let d=0; d < 3; d++){
                        cellPos[d] = i*lattVecs[0][d] + j*lattVecs[1][d] + k*lattVecs[2][d];
                    }

                    for(let a=0; a < atomsPerCell; a++){
                        var atomName = cellTypes[a].replace(/[^a-zA-Z]/gm, "");
                        var t = types.indexOf(cellTypes[a]);
                        var atomPos = cellPos.map((v,d) => v+basisVecs[a][d]);

                        // TODO: Replace each of these repeated geometries with
                        // a single InstancedMesh,
                        // https://threejs.org/docs/#api/en/objects/InstancedMesh.
                        var geometry = new THREE.SphereGeometry(atomRadius, 20, 20);

                        // use closest center for alpha center point
                        var dc = bondCenters.map((v,ii) => mag(diff(atomPos, v)));
                        var c = dc.indexOf(Math.min(...dc));

                        var material = new THREE.MeshPhongMaterial({
                            color: atomColors[atomName],
                            specular: 0x050505,
                            shininess: 100,
                            transparent: true,
                            opacity: alpha(atomPos, bondCenters[c])
                        })

                        // record min/max coordinate for each dimension
                        for(let d=0; d < 3; d++){
                            if(atomPos[d] < coordMinMax[d][0]){
                                coordMinMax[d][0] = atomPos[d];
                            }
                            if(atomPos[d] > coordMinMax[d][1]){
                                coordMinMax[d][1] = atomPos[d];
                            }
                        }

                        atoms[t].push( new THREE.Mesh(geometry, material) );
        
                        scene.add(atoms[t][typeCounts[t]]);
                        atoms[t][typeCounts[t]].position.set(...atomPos);
                        typeCounts[t]++;
                    }
                }
            }
        }
        // coordinate axes
        var O = [0,0,0];
        for(let i=0; i < 3; i++){
            for(let j=0; j < 3; j++){
                O[j] -= lattVecs[i][j]/2;
            } 
        }
        var axesOrigin = new THREE.Vector3(...O);

        var axes = [];
        for(let i=0; i < 3; i++){
            var dir = new THREE.Vector3(...lattVecs[i]);
            dir.normalize();
            var length = mag(lattVecs[i]);

            axes.push( new THREE.ArrowHelper(dir, axesOrigin, length, 0x000000) )
        }

        // trick from Cole's code to generate boundaries for unit cells
        var material = new THREE.LineBasicMaterial({color: "#BEBEBE", opacity: 0.75});
        var lattVecLines = [];
        var pos1 = [0.0, 0.0, 0.0];
        var pos2 = [0.0, 0.0, 0.0];

        for(let j=0; j < lattCells[1]+1; j++){
            for(let k=0; k < lattCells[2]+1; k++){

                for(let d=0; d < 3; d++){
                    pos1[d] = j*lattVecs[1][d] + k*lattVecs[2][d];
                }
                for(let d=0; d < 3; d++){
                    pos2[d] = pos1[d] + lattCells[0]*lattVecs[0][d];
                }

                var geometry = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(...pos1), 
                    new THREE.Vector3(...pos2)
                ]);
                lattVecLines.push(new THREE.Line(geometry, material));
            }
            for(let i=0; i < lattCells[0]+1; i++){

                for(let d=0; d < 3; d++){
                    pos1[d] = j*lattVecs[1][d] + i*lattVecs[0][d];
                }
                for(let d=0; d < 3; d++){
                    pos2[d] = pos1[d] + lattCells[2]*lattVecs[2][d];
                }

                var geometry = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(...pos1), 
                    new THREE.Vector3(...pos2)
                ]);
                lattVecLines.push(new THREE.Line(geometry, material));
            }
        }
        for(let k=0; k < lattCells[2]+1; k++){
            for(let i=0; i < lattCells[0]+1; i++){

                for(let d=0; d < 3; d++){
                    pos1[d] = k*lattVecs[2][d] + i*lattVecs[0][d];
                }
                for(let d=0; d < 3; d++){
                    pos2[d] = pos1[d] + lattCells[1]*lattVecs[1][d];
                }

                var geometry = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(...pos1),
                    new THREE.Vector3(...pos2)
                ]);
                lattVecLines.push(new THREE.Line(geometry, material));
            }
        }

        // make pronounced box around central unit cell where bonds are shown
        // super ugly, but whatever
        material = new THREE.LineBasicMaterial({color: 0x000000, opacity: 1.0});
        var start = centerCell.slice();
        for(let jj=0; jj < 2; jj++){
            var v1 = 0;
            var v2 = 1;
            for(let ii=0; ii < 2; ii++){
                pos1 = start.slice();
                for(let i=0; i < 2; i++){
                    for(let d=0; d < 3; d++){
                        pos2[d] = pos1[d] + lattVecs[v1][d];
                    }
                    var geometry = new THREE.BufferGeometry().setFromPoints([
                        new THREE.Vector3(...pos1),
                        new THREE.Vector3(...pos2)
                    ]);
                    scene.add(new THREE.Line(geometry, material));

                    for(let d=0; d < 3; d++){
                        pos1[d] += lattVecs[v2][d];
                    }
                }
                v1 = 1;
                v2 = 0;
            }
            for(let d=0; d < 3; d++){
                start[d] += lattVecs[2][d];
            }
        }
        start = centerCell.slice();
        for(let jj=0; jj < 2; jj++){
            var v1 = 2;
            var v2 = 1;
            for(let ii=0; ii < 2; ii++){
                pos1 = start.slice();
                for(let i=0; i < 2; i++){
                    for(let d=0; d < 3; d++){
                        pos2[d] = pos1[d] + lattVecs[v1][d];
                    }
                    var geometry = new THREE.BufferGeometry().setFromPoints([
                        new THREE.Vector3(...pos1),
                        new THREE.Vector3(...pos2)
                    ]);
                    scene.add(new THREE.Line(geometry, material));

                    for(let d=0; d < 3; d++){
                        pos1[d] += lattVecs[v2][d];
                    }
                }
                v1 = 1;
                v2 = 2;
            }
            for(let d=0; d < 3; d++){
                start[d] += lattVecs[0][d];
            }
        }

        // convert hex color to RGB string
        function toColor(num) {
            num >>>= 0;
            var b = num & 0xFF,
                g = (num & 0xFF00) >>> 8,
                r = (num & 0xFF0000) >>> 16,
                a = ( (num & 0xFF000000) >>> 24 ) / 255 ;
            return "rgb(" + [r, g, b].join(",") + ")";
        }

        // https://stackoverflow.com/questions/11867545/change-text-color-based-on-brightness-of-the-covered-background-area
        function getContrastYIQ(hexcolor){
            hexcolor = hexcolor.replace("0x", "");
            var r = parseInt(hexcolor.substr(0,2),16);
            var g = parseInt(hexcolor.substr(2,2),16);
            var b = parseInt(hexcolor.substr(4,2),16);
            var yiq = ((r*299)+(g*587)+(b*114))/1000;
            return (yiq >= 128) ? 'black' : 'white';
        }

        // construct bonds represented as tubes
        var bondLines = [];
        for(let i=0; i < bondLabels.length; i++){
            bondLines.push([]);
            for(let j=0; j < bondVecs[i].length; j++){
                var beg = bondCenters[bondTypeIds[i][0]];
                var end = bondVecs[i][j].map((v,k) => v + beg[k]);
                // TODO: Might be cleaner to replace with a CylinderGeometry 
                var geometry = new THREE.TubeGeometry(
                    new THREE.CatmullRomCurve3([ 
                        new THREE.Vector3(...beg),
                        new THREE.Vector3(...end)
                    ]),
                    2, // path segments
                    0.025*maxLattUnit, // thickness
                    12, // roundness
                    false // closed
                );

                const material = new THREE.MeshStandardMaterial({color: toColor(bondColors[i])});
                bondLines[i].push( new THREE.Mesh(geometry, material));
            }
        }

        // create the renderer
        const renderer = new THREE.WebGLRenderer();
        renderer.setSize(1024,1024);
        renderer.setPixelRatio(window.devicePixelRatio); 
        container.append(renderer.domElement);

        function render() {
            renderer.render(scene, camera);
        }

        // Create a camera
        var fov = 0.0 // field of view
        for(let i=0; i < 3; i++){
            fov += (coordMinMax[i][1] - coordMinMax[i][0]);
        }
        const aspect = container.clientWidth / container.clientHeight;
        const near = 0.1; // the near clipping plane
        const far = 20*fov; // far clipping plane
        const camera = new THREE.PerspectiveCamera(fov, aspect, near, far);
        camera.position.set(1.5*coordMinMax[0][1], 1.5*coordMinMax[1][1], 1.5*coordMinMax[2][1]);
        scene.add(camera);

        // create a point light
        var pointLight = new THREE.PointLight( 0xFFFFFF );
        pointLight.position.x = 10;
        pointLight.position.y = 50;
        pointLight.position.z = 130;
        pointLight.intensity = 1.5;
        scene.add(pointLight);
        camera.add(pointLight);

        // orbit controls
        const controls = new OrbitControls( camera, renderer.domElement );
        controls.target.set(...bondCenters[0]);
        controls.addEventListener('change', render);

        // add sublattice options
        var subLattInnerStr = "";
        for(let i=0; i < types.length; i++){
            subLattInnerStr += 
            ` 
            <div>
                <input type="checkbox" id="${types[i]+key}" checked="checked">
                <label for="$types[i] select${key}" style="color: black">${types[i]}</label>
            </div>
            `;
        }
        var subLattDiv = document.createElement("DIV");
        subLattDiv.innerHTML = 
        `<div style="color: black; text-align: left">
            Sublattices
            ${subLattInnerStr}
        </div> &nbsp;
        `;

        document.getElementById(`widget-container-${key}`).appendChild(subLattDiv);
        var subLattToggles = [];
        for(let i=0; i < types.length; i++){
            subLattToggles.push(document.getElementById(types[i]+key));
        }

        // add bond display options
        var bondsInnerStr = "";
        for(let i=0; i < bondLabels.length; i++){
            bondsInnerStr += 
            ` 
            <div style="background-color: ${toColor(bondColors[i])}" id="${bondLabels[i]+key} show">
                <input type="checkbox" id="${bondLabels[i]+key} select" value="0">
                <label for="${bondLabels[i]+key} select" style="color: ${getContrastYIQ(bondColors[i])}">${bondLabels[i]} (${bondTypes[i][0]}, ${bondTypes[i][1]})</label>

                <div style="text-align: left; display: none" id="${bondLabels[i]+key} options"> 
                    <input type="button" id="${bondLabels[i]+key} next" value="${bondVecs[i].length}/${bondVecs[i].length}"/> 
                    <label for="${bondLabels[i]+key} all" style="color: ${getContrastYIQ(bondColors[i])}"><input type="checkbox" id="${bondLabels[i]+key} all">all</label>
                </div>
            </div>
            `;
        }
        var bondsDiv = document.createElement("DIV");
        bondsDiv.innerHTML = 
        ` 
        <div class="row" style=" color: black; text-align: left">
            Bonds
            ${bondsInnerStr}
        </div> &nbsp;
        `;
        document.getElementById(`widget-container-${key}`).appendChild(bondsDiv);

        var selectOne = [];
        var selectAll = [];
        var selectNext = [];
        var bondOptions = [];
        var nextCount = new Array(bondLabels.length).fill(0);

        for(let i=0; i < bondLabels.length; i++){
            selectOne.push( document.getElementById(`${bondLabels[i]+key} select`) );
            bondOptions.push( document.getElementById(`${bondLabels[i]+key} options`) );
            selectAll.push( document.getElementById(`${bondLabels[i]+key} all`) );
            selectNext.push( document.getElementById(`${bondLabels[i]+key} next`) );

            selectAll[i].oninput = function(){
                if(selectAll[i].checked){
                    for(let j=0; j < bondLines[i].length; j++){
                        scene.add(bondLines[i][j]);
                    }
                    selectNext[i].value = `*/${bondVecs[i].length}`;
                }
                else{
                    for(let j=0; j < bondLines[i].length; j++){
                        scene.remove(bondLines[i][j]);
                    }
                    selectNext[i].value = `${nextCount[i]+1}/${bondVecs[i].length}`;
                    scene.add(bondLines[i][nextCount[i]]);
                }
                render()
            }

            selectOne[i].oninput = function(){
                if(selectOne[i].checked){
                    selectAll[i].checked = false;
                    for(let j=0; j < bondLines[i].length; j++){
                        scene.remove(bondLines[i][j]);
                    }
                    scene.add(bondLines[i][nextCount[i]]);
                    bondOptions[i].style.display = "block";
                    selectNext[i].value = `${nextCount[i]+1}/${bondLines[i].length}`;
                }
                else{
                    for(let j=0; j < bondLines[i].length; j++){
                        scene.remove(bondLines[i][j]);
                    }
                    bondOptions[i].style.display = "none";
                }
                render();
            }

            selectNext[i].onclick = function(){
                if(selectAll[i].checked){
                    selectAll[i].checked = false;
                    for(let j=0; j < bondLines[i].length; j++){
                        scene.remove(bondLines[i][j]);
                    }
                }
                scene.remove(bondLines[i][nextCount[i]]);
                nextCount[i] = (++nextCount[i]) % bondLines[i].length;
                scene.add(bondLines[i][nextCount[i]]);
                selectNext[i].value = `${nextCount[i]+1}/${bondLines[i].length}`;
                render();
            }
        }

        // add clear-bonds, unit cell, and axes toggles
        var togglesDiv = document.createElement("DIV");
        togglesDiv.innerHTML = 
        `
        <div>
            <div class="row">  
                <input type="button" id="clear bonds${key}" value="clear bonds"/>
            </div> &nbsp; 

            <div class="row" style="background-color: LightGray; text-align: left">  
                <input type="checkbox" id="axes toggle${key}"/>
                <label for="axes${key} toggle" style="color: black">Show axes</label>
            </div> &nbsp;

            <div class="row" style="background-color: LightGray; text-align: left">  
                <input type="checkbox" id="lattVecs toggle${key}"/>
                <label for="lattVecs${key} toggle" style="color: black">Show unit cells</label>
            </div> &nbsp;
        </div>
        `;
        document.getElementById(`widget-container-${key}`).appendChild(togglesDiv);

        // set clear bonds button functionality
        var clearBonds = document.getElementById("clear bonds"+key);
        clearBonds.onclick = function(){
            for(let i=0; i < bondLines.length; i++){
                for(let j=0; j < bondLines[i].length; j++){
                    scene.remove(bondLines[i][j]);
                }
                selectOne[i].checked = false;
                selectAll[i].checked = false;
                bondOptions[i].style.display = "none";
            }
            render();
        }

        // set lattVecs toggle functionality
        var lattVecsToggle = document.getElementById("lattVecs toggle"+key);

        lattVecsToggle.oninput = function(){
            if(lattVecsToggle.checked){
                for(let j=0; j < lattVecLines.length; j++){
                    scene.add(lattVecLines[j]);
                }
            }
            else{
                for(let j=0; j < lattVecLines.length; j++){
                    scene.remove(lattVecLines[j]);
                }
            }
            render();
        }

        // set axes toggle functionality
        var axesToggle = document.getElementById("axes toggle"+key);

        axesToggle.oninput = function(){
            if(axesToggle.checked){
                for(let j=0; j < axes.length; j++){
                    scene.add(axes[j]);
                }
            }
            else{
                for(let j=0; j < axes.length; j++){
                    scene.remove(axes[j]);
                }
            }
            render();
        }

        // set sublattice toggle functionality
        for(let i=0; i < types.length; i++){
            subLattToggles[i].oninput = function(){
                if(subLattToggles[i].checked){
                    for(let j=0; j < atoms[i].length; j++){
                        scene.add(atoms[i][j]);
                    }

                    for(let j=0; j < bondLabels.length; j++){
                        var flag = true;
                        for(let k=0; k < cellTypes.length; k++){
                            if(bondTypeIds[j].includes(k) && !subLattToggles[types.indexOf(cellTypes[k])].checked){
                                flag = false;
                            }
                        }
                        if(flag){
                            document.getElementById(`${bondLabels[j]} show`+key).style.display = "block";
                        }
                    }
                }
                else{
                    for(let j=0; j < atoms[i].length; j++){
                        scene.remove(atoms[i][j]);
                    }

                    for(let j=0; j < bondLabels.length; j++){
                        var flag = false;
                        for(let k=0; k < bondTypeIds[j].length; k++){
                            if(cellTypes[bondTypeIds[j][k]] == types[i]){
                                flag = true;
                            }
                        }
                        if(flag){
                            selectOne[j].checked = false;
                            selectAll[j].checked = false;
                            for(let k=0; k < bondLines[j].length; k++){
                                scene.remove(bondLines[j][k]);
                            }
                            bondOptions[j].style.display = "none";
                            document.getElementById(`${bondLabels[j]} show`+key).style.display = "none";
                        }
                    }
                }
                render();
            }
        }

        render();
    }

    catch(err) {
        showError(err);
    }
})();
