mapboxgl.accessToken = 'pk.eyJ1Ijoic3NjY2hhbiIsImEiOiJjajB6bm5zOTkwMnI3MndwbXhieXZ2a3c3In0.Zg6_NdjbijgrnNLfT-SJEQ';

var tileset = 'mapbox.streets';
var map = new mapboxgl.Map({
    container: 'map', // container id
    style: {
        "version": 8,
        "sources": {
            "raster-tiles": {
                "type": "raster",
                "url": "mapbox://" + tileset,
                "tileSize": 256
            }
        },
        "layers": [{
            "id": "simple-tiles",
            "type": "raster",
            "source": "raster-tiles",
            "minzoom": 0,
            "maxzoom": 22
        }]
    },
    center: [46.460938, -19.002846], // starting position
    zoom: 5, // starting zoom
});

map.on('style.load', function() {
  addAllDataLayers();
  renderLayers();
});

var visibility = {
   "radiance-layer": false,
   "population-layer": false,
   "viableElectrification-layer": false,
   "transmission-layer": false,
   "road-layer": false,
};

function renderLayers() {
  for (layer in visibility) {
    if (visibility[layer] === true) {
      map.setLayoutProperty(layer, 'visibility', 'visible');
    } else {
      map.setLayoutProperty(layer, 'visibility', 'none');
    }
  }
}

function addAllDataLayers() {
  addLightPolutionLayer();
  addPopulationHeatmapLayer();
  addViableElectrification();  
  addTransmissionLayer();
  addRoadLayer();
}

// Data Sources
function addLightPolutionLayer() {
    map.addSource('radiance', {
        type: 'raster',
        url: 'mapbox://sscchan.crg6g9hu'
    });
    map.addLayer({
            "id": "radiance-layer",
            "type": "raster",
            "source": "radiance",
            "minzoom": 0,
            "maxzoom": 22
    });
  
  map.setLayoutProperty('radiance-layer', 'visibility', 'none');
}

function addPopulationHeatmapLayer() {
    map.addSource('population-heatmap', {
        type: 'raster',
        url: 'mapbox://sscchan.0bnrmt0b'
    });
  
    map.addLayer({
            "id": "population-layer",
            "type": "raster",
            "source": "population-heatmap",
            "minzoom": 0,
            "maxzoom": 22
    });

    map.setPaintProperty('population-layer', 'raster-opacity', 0.5);
    map.setLayoutProperty('population-layer', 'visibility', 'none');  
}

function addTransmissionLayer() {
    map.addLayer({
        "id": "transmission-layer",
        "type": "line",
        "source": {
            type: 'vector',
            url: 'mapbox://sscchan.799cga7b'
        },
        "source-layer": "madagascarTransmission2006-d4acqr",
        "layout": {
            "line-join": "round",
            "line-cap": "round"
        },
        "paint": {
            "line-color": "#ff0000",
            "line-width": 2
        },
    });

    map.setLayoutProperty('transmission-layer', 'visibility', 'none');    
}

function addViableElectrification() {
    map.addLayer({
        "id": "viableElectrification-layer",
        "type": "line",
        "source": {
            type: 'vector',
            url: 'mapbox://sscchan.56ez7nj4'
        },
        "source-layer": "settlement_polygon-ctfukj",
        "layout": {
            "line-join": "round",
            "line-cap": "round"
        },
        "paint": {
            "line-color": "#555555",
            "line-width": 2
        },
    });
  
    map.setLayoutProperty('road-layer', 'visibility', 'none');
}

function addRoadLayer() {
    map.addLayer({
        "id": "road-layer",
        "type": "line",
        "source": {
            type: 'vector',
            url: 'mapbox://sscchan.bmwrmc23'
        },
        "source-layer": "madagascar-roads-82554w",
        "layout": {
            "line-join": "round",
            "line-cap": "round"
        },
        "paint": {
            "line-color": "#000000",
            "line-width": 2
        },
    });
  
    map.setLayoutProperty('road-layer', 'visibility', 'none');
}


// Event Handlers
$(document).ready(function() {
  $("#radianceButton").on("click", function(){
    visibility["radiance-layer"] = !visibility["radiance-layer"];
    renderLayers();
  });
  $("#populationButton").on("click", function(){
    visibility["population-layer"] = !visibility["population-layer"];
    renderLayers();
  });
  $("#transmissionButton").on("click", function(){
    visibility["transmission-layer"] = !visibility["transmission-layer"];
    renderLayers();
  });  
  $("#roadsButton").on("click", function(){
    visibility["road-layer"] = !visibility["road-layer"];
    renderLayers();
  });    

  $("#viableElectrificationButton").on("click", function(){
    visibility["viableElectrification-layer"] = !visibility["viableElectrification-layer"];
    renderLayers();
  });    
  
})

// Search Box
map.addControl(new MapboxGeocoder({
    accessToken: mapboxgl.accessToken
}));

map.addControl(new mapboxgl.ScaleControl({
    maxWidth: 250,
    unit: 'metric'
}));

// Panel to toggle satellite view and other views

var layerList = document.getElementById('menu');
var inputs = layerList.getElementsByTagName('input');

function switchLayer(layer) {
    var layerId = layer.target.id;
    map.setStyle('mapbox://styles/mapbox/' + layerId + '-v9');
}

for (var i = 0; i < inputs.length; i++) {
    inputs[i].onclick = switchLayer;
}