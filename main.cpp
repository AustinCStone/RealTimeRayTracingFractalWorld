#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "gl.h"

Shader shader;
Buffer<vec2> quad;
VAO layout;

void setup() {
    shader.vertexShader(glsl(
      attribute vec2 vertex;
      varying vec2 screenPosition;

      void main() {
        screenPosition = (vertex+.5);
        gl_Position = vec4(vertex, 0.0, 1.0);
      }
   )).fragmentShader(glsl(
      uniform vec3 eye;
      uniform vec3 axisX;
      uniform vec3 axisY;
      uniform vec3 axisZ;
      uniform float angleOfRotation;
      uniform vec3 lightDirection;
      uniform float maxVisibilityDueToFog;
      uniform int enableMoon;
      uniform int enableRipples;
      uniform int enableSnow;
      uniform int enableClouds;
      uniform float timeElapsed;
      uniform float currentJumpHeight;


      varying vec2 screenPosition;

      const int NUM_GRIDS = 5;
      float GRID_WEIGHTINGS[NUM_GRIDS];
      float GRID_SPACINGS[NUM_GRIDS];

      // Material Properties
      const int NUM_MATERIALS = 3;
      const int MOUNTAIN = 0;
      const int WATER = 1;
      const int SNOW = 2;

      float diffuseCoefficients[NUM_MATERIALS];
      float specularCoefficients[NUM_MATERIALS];
      float specularExponents[NUM_MATERIALS];
      float reflectivityCoefficients[NUM_MATERIALS];
      bool hasReflection[NUM_MATERIALS];

      //vec3  LIGHT_DIRECTION = normalize(vec3(1.0, 0.5, 1.0));
      int ALIASING_AMOUNT = 1;
      float ALIASING_EPSILON = .00012;
      float CAMERA_ELEVATION = 5.0;
      float CAM_DISTANCE_TO_SCREEN = 0.05;
      float WORLD_VIEWING_WIDTH = 0.05;
      float RAY_EPSILON = 0.001;
      float MAX_HEIGHT = 25.0;
      float ASPECT_RATIO = 0.625;
      float NEIGHBOR_DISTANCE = 0.01;
      float DIFFUSE_COEF = 0.7;
      float SPECULAR_COEF = 0.3;
      float SPECULAR_EXPONENT = 1.0;
      float LIGHT_BRIGHTNESS = 1.0;
      //float MAX_HORIZON = 50.0;
      // float REFLECTION_COEF = 0.1;
      vec3 FOG_COLOR = vec3(.2, .4, .5);
      float SUN_SIZE = 0.13;
      float MIN_INTENSITY_OF_SUN_APPEARANCE = 0.3;  
      float WATER_HEIGHT = 10.0;
      float WATER_SPECULAR_EXPONENT = 3.0;
      float CLOUD_CEIL = 18.0;
      float CLOUD_FLOOR = 15.0;

	    // Set via the keyboard ("m" key)
      bool ENABLE_MOON;
      // Set via the keyboard ("s" key)
      bool ENABLE_SNOW;
      // Set via the keyboard ("h" key)
      bool ENABLE_RIPPLES;
      // Set via the keyboard ("c" key)
      bool ENABLE_CLOUDS;
      

      float smoothStep(in float x)
      {
        return x * x * (3.0 - (2.0 * x));
      }

      //borrowed from http://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
      float hash( float n )
      {
        return fract(sin(n)*43758.5453);
      }


      //borrowed from http://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
      float cloudNoise( vec3 x )
      {
        // The noise function returns a value in the range -1.0f -> 1.0f

        vec3 p = floor(x);
        vec3 f = fract(x);

        f = f*f*(3.0-2.0*f);
        float n = p.x + p.y*57.0 + 113.0*p.z;

        return pow((mix(mix(mix( hash(n+0.0), hash(n+1.0),f.x),
                   mix( hash(n+57.0), hash(n+58.0),f.x),f.y),
               mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                   mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z)) / 0.5, 4.0) / 10.0;
      }


      float attenuateCloudDensitySample(in vec3 location, in float density)
      {
        float distanceFromCenter = abs(location.y - ((CLOUD_FLOOR + CLOUD_CEIL)/2.0));
        return density * smoothStep((1.0 - (distanceFromCenter / ((CLOUD_CEIL - CLOUD_FLOOR)/2.0))));
      }

      float getCloudDensitySample(in vec3 location) 
      {
        float density = cloudNoise( location - vec3(0.0,0.1,1.0) * 20.0 * timeElapsed ); 
        return attenuateCloudDensitySample(location, density);
      }

      vec4 integrateCloudRay(in vec4 currentCloudColor, in float nextSampleDensity)
      {

        vec4 nextSampleColor = vec4( mix(vec3(1.0, .4 ,0.4), vec3(0.4, .5,0.9), nextSampleDensity*2.5 ), nextSampleDensity );
        nextSampleColor.rgb *= nextSampleColor.a;
        return currentCloudColor + nextSampleColor*(1.0-currentCloudColor.a);
      }

      float getDistance(in vec3 start, in vec3 end) 
      {
        return sqrt(dot((end-start), (end-start)));
      }

      vec4 getCloudContribution(in vec3 rayStartingLocation, in vec3 rayDirection, in vec3 rayEndLocation)
      {
        vec3 viewerPosition = rayStartingLocation;

        if(rayStartingLocation.y < CLOUD_FLOOR)
        {
          rayStartingLocation = rayStartingLocation + rayDirection * ((CLOUD_FLOOR - rayStartingLocation.y) / rayDirection.y);
        }

        if(rayStartingLocation.y > CLOUD_CEIL)
        {

          rayStartingLocation = rayStartingLocation + rayDirection * ((CLOUD_CEIL - rayStartingLocation.y) / rayDirection.y);
        }

        float distanceToStartOfClouds = getDistance(viewerPosition, rayStartingLocation);
        float attenuationOfClouds = min(max(1.0 - (distanceToStartOfClouds / maxVisibilityDueToFog), 0.0), 1.0);

        if(rayEndLocation.y < CLOUD_FLOOR && rayEndLocation.y > -100.0)
        {
          rayEndLocation = rayEndLocation + rayDirection * ((CLOUD_FLOOR - rayEndLocation.y) / rayDirection.y);
        }

        if(rayEndLocation.y > CLOUD_CEIL || rayEndLocation.y < -100.0)
        {
          if (rayEndLocation.y<-100.0) {

           
            if(rayDirection.y< 0.0) {
              rayEndLocation = rayStartingLocation + rayDirection * ((CLOUD_FLOOR - rayStartingLocation.y) / rayDirection.y);
            } else {
              rayEndLocation = rayStartingLocation + rayDirection * ((CLOUD_CEIL - rayStartingLocation.y) / rayDirection.y);
            }
          } else {
              
            rayEndLocation = rayEndLocation + rayDirection * ((CLOUD_CEIL - rayEndLocation.y) / rayDirection.y);
          }
        }






        //March from rayStartingLocation to rayEndLocation, and return the cloud color and density of the ray  
        vec3 currentRayPosition = rayStartingLocation;
        vec4 cloudSum = vec4(0.0, 0.0, 0.0, 0.0);
        //float t_increment = 0.1;
        float t = 0.0;
        float endT = (rayEndLocation.x-rayStartingLocation.x)/rayDirection.x;
        while(t<endT && cloudSum.a < .99)
        {
          float nextSampleDensity = getCloudDensitySample(currentRayPosition);
          cloudSum = integrateCloudRay(cloudSum, nextSampleDensity);
          t += max(0.05,0.02*t); 
          currentRayPosition = rayStartingLocation + (t * rayDirection);
        }
       /* float rayDistanceInSector = getDistance(rayStartingLocation, rayEndLocation);
        if(rayDistanceInSector>(CLOUD_CEIL-CLOUD_FLOOR)){
          return cloudSum * (CLOUD_CEIL-CLOUD_FLOOR)/rayDistanceInSector;
        }
        return cloudSum;*/
        return cloudSum * attenuationOfClouds; //* abs(dot(rayDirection, vec3(0.0, 1.0, 0.0)));
      }

      //Hash12 Borrowed from David Hospkins Shader Toy Demo Code
      float Hash12(in float x, in float y)
      {
        float xSave = 3.07965;
        float ySave = 7.4235;

        x = (x / xSave) - floor(x / xSave);
        y = (y / ySave) - floor(y / ySave);

        float xUpdate = y + (x * y) + (y * x) + 19.19;
        float yUpdate = x + (x * y) + (y * x) + 19.19;

  		  float value = MAX_HEIGHT * ((xUpdate * yUpdate) - floor(xUpdate * yUpdate));
  		  
  		  if(ENABLE_MOON)
  		  {
  			  // Moon-Like Surface Code  
  			  // Make valleys, not mountains
  			  value = -1.0 * value;
  			  if(abs(value) < 22.0)
  			  {
  			  	// Level out parts of terrain that are relatively close to y = 0 to begin with
  				  value = value * 0.25;
  			  }
  			  else
  			  {
  			  	// Make the valleys not quite as sizeable as the mountains would be
  			  	value = value * 0.75;
  			  }
  		  }
            return value; 
      }


      float getFogginess(in vec3 intersectionLocation, in vec3 rayStartingLocation) {
        return  max(min(sqrt(dot(intersectionLocation-rayStartingLocation, intersectionLocation-rayStartingLocation))/maxVisibilityDueToFog, 1.0),0.0);
      }

      



      float interpolateSimple(in float xGridFract, in float zGridFract, in float xGridCoord, in float zGridCoord)
      {
       
        xGridFract = smoothStep(xGridFract);
        zGridFract = smoothStep(zGridFract);

        float topLeft = Hash12(xGridCoord, zGridCoord);
        float topRight = Hash12(xGridCoord + 1.0, zGridCoord);
        float bottomLeft = Hash12(xGridCoord, zGridCoord + 1.0);
        float bottomRight = Hash12(xGridCoord + 1.0, zGridCoord + 1.0);
    

        float topLeftContribution = (1.0 - xGridFract) * (1.0 - zGridFract) * topLeft;
        float topRightContribution = xGridFract * (1.0 - zGridFract) * topRight;
        float bottomLeftContribution = (1.0 - xGridFract) * zGridFract * bottomLeft;
        float bottomRightContribution = xGridFract * zGridFract * bottomRight;
 
        return topLeftContribution + topRightContribution + bottomLeftContribution + bottomRightContribution;
      }


      float getHeightByLocation(in float x, in float z)
      {
     
        float height = 0.0;

        int grid_element = 0;

        for(int i = 0; i < NUM_GRIDS; i++)
        {
          float gridSpacing = GRID_SPACINGS[i];
          float gridWeighting = GRID_WEIGHTINGS[i];

          float xGridCoord = floor(x / gridSpacing);
          float zGridCoord = floor(z / gridSpacing);

          float xGridFract = (x - (xGridCoord * gridSpacing)) / gridSpacing;
          float zGridFract = (z - (zGridCoord * gridSpacing)) / gridSpacing;

          height += (interpolateSimple(xGridFract, zGridFract, xGridCoord, zGridCoord) * gridWeighting);
        }
        return height;
      }


      vec3 advanceRay(in vec3 rayDirection, in vec3 rayStartingLocation, float tDelta)
      {
        return rayStartingLocation + (rayDirection * tDelta);
      }


      vec3 binarySearchIntersection(in float startT, in float endT, in vec3 rayDirection, in vec3 rayStartingLocation)
      {
        vec3 currentRayLocation = advanceRay(rayDirection, rayStartingLocation, startT);
        while((endT - startT) > .00001)
        {
          float midT = (startT / 2.0) + (endT / 2.0);
          currentRayLocation = advanceRay(rayDirection, rayStartingLocation, midT);
          float heightDelta = currentRayLocation.y - getHeightByLocation(currentRayLocation.x, currentRayLocation.z);

          if (heightDelta > 0.0)
          {
            startT = midT;
          }
          else
          {
            endT = midT;
          }
        }
        return currentRayLocation;
      }


      vec3 getRayIntersection(in vec3 rayDirection, in vec3 rayStartingLocation)
      {
          float tDelta = 0.05;
          vec3 currentRayLocation = rayStartingLocation;
          float t = 0.0;

          float heightDelta = currentRayLocation.y - getHeightByLocation(currentRayLocation.x, currentRayLocation.z);

          while(heightDelta > RAY_EPSILON)
          {
            tDelta = .3 * heightDelta;
            t += tDelta;

            if (t > 50.0)
            {
              // TODO: Change to a different flag?
              return vec3(0.0, -1000.0, 0.0);
            }

            currentRayLocation = advanceRay(rayDirection, currentRayLocation, tDelta);
            heightDelta = currentRayLocation.y - getHeightByLocation(currentRayLocation.x, currentRayLocation.z);
          }

          return binarySearchIntersection(t - (2.0*tDelta), t + tDelta, rayDirection, rayStartingLocation);
      }


      vec3 getNormalAtLocation(in vec3 location, in float scalingOfNeighborDistance)
      {
        float neighborHeightToRight = getHeightByLocation(location.x + NEIGHBOR_DISTANCE * scalingOfNeighborDistance, location.z);
        float neighborHeightInFront = getHeightByLocation(location.x, location.z + NEIGHBOR_DISTANCE * scalingOfNeighborDistance);

        vec3 positionOfRightNeighbor = vec3(location.x + NEIGHBOR_DISTANCE * scalingOfNeighborDistance, neighborHeightToRight, location.z);
        vec3 positionOfFrontNeighbor = vec3(location.x, neighborHeightInFront, location.z + NEIGHBOR_DISTANCE * scalingOfNeighborDistance);

        vec3 directionToRightNeighbor = normalize(positionOfRightNeighbor - location);
        vec3 directionToFrontNeighbor = normalize(positionOfFrontNeighbor - location);

        vec3 normalAtCurrentLocation = normalize(cross(directionToFrontNeighbor, directionToRightNeighbor));
        return normalAtCurrentLocation;
      }


      float getSpecularContribution(in vec3 locationInWorld, in vec3 normalAtLocation, in vec3 viewingDirection, in int terrainMaterial)
      {
        vec3 incomingLightRay = normalize( -1.0 * lightDirection);
        vec3 reflectionRay = normalize(reflect(incomingLightRay, normalAtLocation));

        float dotProduct = max(dot(viewingDirection, reflectionRay),0.0);

        float specularContribution = pow(dotProduct, specularExponents[terrainMaterial]);

        return LIGHT_BRIGHTNESS * specularCoefficients[terrainMaterial] * specularContribution;
      }


      float getDiffuseContribution(in vec3 normalAtLocation, in int terrainMaterial) 
      {
        return LIGHT_BRIGHTNESS * max(diffuseCoefficients[terrainMaterial] * dot(normalize(normalAtLocation), normalize(lightDirection)), 0.0);
      }


      bool hasRayToLight(in vec3 location) 
      {
        vec3 intersectionLocation = getRayIntersection(lightDirection, location + (lightDirection * 4.0 * RAY_EPSILON));
        if(intersectionLocation.y < -100.0)
        {
          return true;
        }
        else
        {
          return false;
        }
      }
      

      float calculateAngleBetweenRayAndSun(in vec3 ray)
      {
      	float dotProduct = dot(ray, lightDirection); 	
      	float productOfLengths = length(ray) * length(lightDirection);
      	float cosineOfTheAngle = dotProduct / productOfLengths;
      	return acos(cosineOfTheAngle);
      }
      

      vec3 getColorOfSun(in float angleFromCenterOfSun)
      {
      	float intensity = 1.0 - (angleFromCenterOfSun / SUN_SIZE);
        intensity = smoothStep(intensity);
            	
        // The sun can be up to twice as bright as its minimum brightness, depending on the amount of fog
        float intensityAmount = MIN_INTENSITY_OF_SUN_APPEARANCE + MIN_INTENSITY_OF_SUN_APPEARANCE * (maxVisibilityDueToFog / 50.0);
            
        return (intensity * vec3(FOG_COLOR.x + intensityAmount, FOG_COLOR.y + intensityAmount, FOG_COLOR.z + intensityAmount)) + (1.0 - intensity) * FOG_COLOR;
      }
      

      vec3 getWaterIntersectionLocation(in vec3 rayStartingLocation, in vec3 rayDirection)
      {
        float scaleFactorToWater = (WATER_HEIGHT - rayStartingLocation.y) / rayDirection.y;
        return vec3(rayStartingLocation + rayDirection * scaleFactorToWater);
      }
      

      // http://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
      float rand(vec2 co)
      {
        return fract(sin(dot(co.xy ,co.xy)));
   	  }


      vec3 getReflectionContribution(in vec3 rayStartingLocation, in vec3 rayDirection) 
      {
        float totalLightRed = 0.0;
        float totalLightGreen = 0.0;
        float totalLightBlue = 0.0;
        
        vec3 intersectionLocation = getRayIntersection(rayDirection, rayStartingLocation);


        vec4 cloudContribution = vec4(0.0, 0.0, 0.0, 0.0);
        //anyttime you are inside the cloud realm

        if(!ENABLE_MOON && ENABLE_CLOUDS)
        {
          if(rayStartingLocation.y>CLOUD_FLOOR && rayStartingLocation.y<CLOUD_CEIL ) {
             cloudContribution = getCloudContribution(rayStartingLocation, rayDirection, intersectionLocation);
          }
          else if(rayStartingLocation.y < CLOUD_FLOOR && (intersectionLocation.y > CLOUD_FLOOR || intersectionLocation.y < -100.0))
          {
              cloudContribution = getCloudContribution(rayStartingLocation, rayDirection, intersectionLocation);
          }
          else if(rayStartingLocation.y > CLOUD_CEIL && intersectionLocation.y < CLOUD_CEIL)
          {
              cloudContribution = getCloudContribution(rayStartingLocation, rayDirection, intersectionLocation);
          }
        }
        
        vec3 normalAtCurrentLocation = getNormalAtLocation(intersectionLocation, 1.0);
        int material = MOUNTAIN;

        //TODO put these values in globals 
        if(ENABLE_SNOW)
        {
          if(intersectionLocation.y > MAX_HEIGHT*.6 && ((normalAtCurrentLocation.y*MAX_HEIGHT/2.5) + intersectionLocation.y/1.0) > MAX_HEIGHT/1.1)
          {
            material = SNOW;
            normalAtCurrentLocation = getNormalAtLocation(intersectionLocation, 12.0);
          }
        }
       

        float specularContriution = getSpecularContribution(intersectionLocation, normalAtCurrentLocation, rayDirection*-1.0, material);
        

        vec3 colorOfReflection = vec3(0.0, 0.0, 0.0);
      
        float fogAttenuation = getFogginess(intersectionLocation, rayStartingLocation);
        float shadowAttenuation = 1.0;
        
        if(! hasRayToLight(intersectionLocation)) 
        {
              shadowAttenuation = 0.5; //make colors half as bright if there is a shadow
        }
        if(intersectionLocation.y < -100.0)
        {   
          float angleBetweenRayAndSun = calculateAngleBetweenRayAndSun(rayDirection);
          if(angleBetweenRayAndSun < SUN_SIZE)
          {
            return getColorOfSun(angleBetweenRayAndSun);
          }
          else
          {
            totalLightRed  += 1.0 * shadowAttenuation; //+ fogAttenuation;
            totalLightGreen += 0.0 * shadowAttenuation;// + fogAttenuation;
            totalLightBlue +=  0.0 * shadowAttenuation; //+ fogAttenuation;
          }
        }
        else
        {
          //Leave this line here for it to run on 2010 mac book pro? 
          float diffuseContribution = getDiffuseContribution(normalAtCurrentLocation, material);
          totalLightRed  += (diffuseContribution + specularContriution) * shadowAttenuation; 
          totalLightGreen += (diffuseContribution + specularContriution) * shadowAttenuation;
          totalLightBlue += (diffuseContribution + specularContriution) * shadowAttenuation;

          if(material == SNOW)
          {
            float snowBrightness = 0.2;
            totalLightRed += snowBrightness;
            totalLightGreen += snowBrightness;
            totalLightBlue += snowBrightness;
          }
        }
          return (1.0 - cloudContribution.a) * (fogAttenuation * vec3(FOG_COLOR.xyz) + (1.0-fogAttenuation) * vec3(totalLightRed, totalLightGreen, totalLightBlue)) + cloudContribution.rgb;

      }

      vec3 getWaterContribution(in vec3 rayStartingLocation, in vec3 rayDirection, in vec3 intersectionLocation)
      {
      	vec3 waterColor = vec3(0.0, 0.0, 0.0);
      	
      	vec3 waterIntersectionLocation = getWaterIntersectionLocation(rayStartingLocation, rayDirection);
        float xValueOfNorm = 0.0;
        float zValueOfNorm = 0.0; 
      	if(ENABLE_RIPPLES) 
        {
          xValueOfNorm = rand(vec2(waterIntersectionLocation.x+timeElapsed, waterIntersectionLocation.z+timeElapsed));
          zValueOfNorm = rand(vec2(waterIntersectionLocation.x+2.0*timeElapsed, waterIntersectionLocation.z+2.0*timeElapsed));
        }
      

		    // The normal sticks up out of the water, but is perturbed to slightly to get a non-uniform surface
        vec3 normalAtWater = normalize(vec3(xValueOfNorm / 120.0, 1.0, zValueOfNorm / 120.0));

        float specularContribution = getSpecularContribution(waterIntersectionLocation, normalAtWater, rayDirection*-1.0, WATER);
        float diffuseContribution = getDiffuseContribution(normalAtWater, WATER);
        float fogAttenuation = getFogginess(waterIntersectionLocation, rayStartingLocation);

        float shadowAttenuation = 1.0;
        
        if(! hasRayToLight(waterIntersectionLocation+vec3(xValueOfNorm/4.0, 0.0, zValueOfNorm/4.0))) 
        {
              shadowAttenuation = 0.85;
        }

        vec3 reflectionContribution = vec3(0.0, 0.0, 0.0);
        if( hasReflection[WATER]) 
        {
          reflectionContribution = getReflectionContribution(waterIntersectionLocation, normalize(reflect(rayDirection, normalAtWater)));
        }
      	
      	waterColor[0] += (specularContribution +   (FOG_COLOR.r *  diffuseContribution) + reflectivityCoefficients[WATER] * reflectionContribution.r) * shadowAttenuation;
        waterColor[1] += (specularContribution +   (FOG_COLOR.g *  diffuseContribution) +  reflectivityCoefficients[WATER] * reflectionContribution.g) * shadowAttenuation;
        waterColor[2] += (specularContribution +  (FOG_COLOR.b *  diffuseContribution) + reflectivityCoefficients[WATER] * reflectionContribution.b) * shadowAttenuation;

        return fogAttenuation * FOG_COLOR + (1.0-fogAttenuation) * waterColor; 
      }

      
  	  vec3 getRayColor(in vec3 rayStartingLocation, in vec3 rayDirection)
      {
          float totalLightRed = 0.0;
          float totalLightGreen = 0.0;
          float totalLightBlue = 0.0;
          
          vec3 intersectionLocation = getRayIntersection(rayDirection, rayStartingLocation);
          vec4 cloudContribution = vec4(0.0, 0.0, 0.0, 0.0);
          //anyttime you are inside the cloud realm

          if(!ENABLE_MOON && ENABLE_CLOUDS)
          {
            if(rayStartingLocation.y>CLOUD_FLOOR && rayStartingLocation.y<CLOUD_CEIL ) {
                cloudContribution = getCloudContribution(rayStartingLocation, rayDirection, intersectionLocation);
            }
            else if(rayStartingLocation.y < CLOUD_FLOOR && (intersectionLocation.y > CLOUD_FLOOR || intersectionLocation.y < -100.0))
            {
              cloudContribution = getCloudContribution(rayStartingLocation, rayDirection, intersectionLocation);
            }
            else if(rayStartingLocation.y > CLOUD_CEIL && intersectionLocation.y < CLOUD_CEIL)
            {
              cloudContribution = getCloudContribution(rayStartingLocation, rayDirection, intersectionLocation);
            }
          }

          if(!ENABLE_MOON)
          {
  		      vec3 waterContribution = vec3(0.0, 0.0, 0.0);
            if(intersectionLocation.y < WATER_HEIGHT && intersectionLocation.y > -900.0)
            {
              waterContribution = getWaterContribution(rayStartingLocation, rayDirection, intersectionLocation);
              return (1.0 - cloudContribution.a) * waterContribution + cloudContribution.rgb;
            }
  		    }
          
          vec3 normalAtCurrentLocation = getNormalAtLocation(intersectionLocation, 1.0);
          int material = MOUNTAIN;

          //TODO put these values in globals 
          if(ENABLE_SNOW)
          {
            if(intersectionLocation.y > MAX_HEIGHT*.6 && ((normalAtCurrentLocation.y*MAX_HEIGHT/2.5) + intersectionLocation.y/1.0) > MAX_HEIGHT/1.1)
            {
              material = SNOW;
              normalAtCurrentLocation = getNormalAtLocation(intersectionLocation, 12.0);
            }
          }
         

          float specularContribution = getSpecularContribution(intersectionLocation, normalAtCurrentLocation, rayDirection*-1.0, material);
          float diffuseContribution = getDiffuseContribution(normalAtCurrentLocation, material);

          vec3 colorOfReflection = vec3(0.0, 0.0, 0.0);
          if(!ENABLE_MOON)
          {
            if(hasReflection[material] && !(intersectionLocation.y < -100.0))
            {
              colorOfReflection = getReflectionContribution(intersectionLocation, normalize(reflect(rayDirection, normalAtCurrentLocation)));
            }
          }
        
          float fogAttenuation = getFogginess(intersectionLocation, rayStartingLocation);
  		    float shadowAttenuation = 1.0;
          
          if(! hasRayToLight(intersectionLocation)) 
          {
                shadowAttenuation = 0.5; //make colors half as bright if there is a shadow
          }
          if(intersectionLocation.y < -100.0)
          {   
            float angleBetweenRayAndSun = calculateAngleBetweenRayAndSun(rayDirection);
            if(angleBetweenRayAndSun < SUN_SIZE)
            {
  			      return  (1.0 - cloudContribution.a) * getColorOfSun(angleBetweenRayAndSun) + cloudContribution.rgb;
            }
            else
            {
            	totalLightRed  += 1.0 * shadowAttenuation; //+ fogAttenuation;
            	totalLightGreen += 0.0 * shadowAttenuation;// + fogAttenuation;
            	totalLightBlue +=  0.0 * shadowAttenuation; //+ fogAttenuation;
            }
          }
          else
          {
            totalLightRed  += ((1.0 - reflectivityCoefficients[material]) *(diffuseContribution + specularContribution) + (reflectivityCoefficients[material] * colorOfReflection[0])) * shadowAttenuation; //fogAttenuation;
            totalLightGreen += ((1.0 - reflectivityCoefficients[material]) * (diffuseContribution + specularContribution) + (reflectivityCoefficients[material] * colorOfReflection[1])) * shadowAttenuation;//fogAttenuation;
            totalLightBlue += ((1.0 - reflectivityCoefficients[material]) * (diffuseContribution + specularContribution) + (reflectivityCoefficients[material] * colorOfReflection[2])) * shadowAttenuation;// fogAttenuation;

            

            if(material == SNOW)
            {
              float snowBrightness = 0.2;
              totalLightRed += snowBrightness;
              totalLightGreen += snowBrightness;
              totalLightBlue += snowBrightness;
            }


          }
          return (1.0 - cloudContribution.a) * (fogAttenuation * vec3(FOG_COLOR.xyz) + (1.0-fogAttenuation) * vec3(totalLightRed, totalLightGreen, totalLightBlue)) + cloudContribution.rgb;
        }
        

        void initializeGridWeightings()
        { 
        	if(ENABLE_MOON)
        	{
        		GRID_WEIGHTINGS[0] = 0.3;
      			GRID_WEIGHTINGS[1] = 0.5;
      			GRID_WEIGHTINGS[2] = 0.1; 
      			GRID_WEIGHTINGS[3] = .003; 
      			GRID_WEIGHTINGS[4] = .001;
      			GRID_SPACINGS[0] = 15.0;
      			GRID_SPACINGS[1] = 9.0;
      			GRID_SPACINGS[2] = 1.0;
      			GRID_SPACINGS[3] = .1;
      			GRID_SPACINGS[4] = .03;
        	}
        	else
        	{
      			GRID_WEIGHTINGS[0] = 0.5;
      			GRID_WEIGHTINGS[1] = 0.5;
      			GRID_WEIGHTINGS[2] = 0.1; 
      			GRID_WEIGHTINGS[3] = .003; 
      			GRID_WEIGHTINGS[4] = .001;
      			GRID_SPACINGS[0] = 20.0;
      			GRID_SPACINGS[1] = 7.0;
      			GRID_SPACINGS[2] = 1.0;
      			GRID_SPACINGS[3] = .1;
      			GRID_SPACINGS[4] = .03;
        	}
        }
        
        void initializeTerrainParameters()
        {
          diffuseCoefficients[MOUNTAIN] = 0.7;
          diffuseCoefficients[WATER] = 0.7;
          diffuseCoefficients[SNOW] = 3.0;

          specularCoefficients[MOUNTAIN] = 0.3;
          specularCoefficients[WATER] = 0.3;
          specularCoefficients[SNOW] = 0.0;

          specularExponents[MOUNTAIN] = 40.0;
          specularExponents[WATER] = 4.0;
          specularExponents[SNOW] = 10.0;

          reflectivityCoefficients[MOUNTAIN] = 0.0;
          reflectivityCoefficients[WATER] = 0.8;
          reflectivityCoefficients[SNOW] = .5;

          hasReflection[MOUNTAIN] = false;
          hasReflection[WATER] = true;
          hasReflection[SNOW] = false;


        	if(enableMoon == 1)
        	{
        		ENABLE_MOON = true;
        	}
        	else
        	{
        		ENABLE_MOON = false;
        	}

          if(enableRipples == 1) 
          {
            ENABLE_RIPPLES = true;
          }
          else
          {
            ENABLE_RIPPLES = false;
          }

          if(enableSnow == 1) 
          {
            ENABLE_SNOW = true;
          }
          else
          {
            ENABLE_SNOW = false;
          }

          if(enableClouds == 1) 
          {
            ENABLE_CLOUDS = true;
          }
          else
          {
            ENABLE_CLOUDS = false;
          }
        
        	if(ENABLE_MOON)
        	{
        		CAMERA_ELEVATION = 2.0;
        		FOG_COLOR = vec3(0.5, 0.5, 0.5);
        	}
        }
        
        float getCameraHeight()
        {
  	  	  if(ENABLE_MOON)
  	  	  {
  	  		 return CAMERA_ELEVATION + getHeightByLocation(eye.x, eye.z) + currentJumpHeight;
  	  	  }
  	  	  else
  	  	  {
            return max(WATER_HEIGHT + 3.0, CAMERA_ELEVATION + getHeightByLocation(eye.x, eye.z) + currentJumpHeight);
  	  	  }
        }

        void main() 
        {
          initializeTerrainParameters();
        
  	  	  initializeGridWeightings();
  	  	
        	vec3 focalLocation = vec3(eye.x, getCameraHeight(), eye.z);

        	float totalLightRed =  0.0;
        	float totalLightGreen =  0.0;
        	float totalLightBlue = 0.0;

        	mat2 rotationMat = mat2(
            cos(angleOfRotation), sin(angleOfRotation),// first column (not row!)
            -sin(angleOfRotation), cos(angleOfRotation)); // second column
          
  	      float fogAttenuation = 0.0;

          for(int i = 0; i<ALIASING_AMOUNT; i++) 
          {
            for(int j = 0; j<ALIASING_AMOUNT; j++) 
            {
              vec3 rayDirection = vec3(0.0, 0.0, 0.0);
      				vec3 rayStartingLocation = vec3(0.0, 0.0, 0.0);

      				rayStartingLocation.x = (screenPosition.x * WORLD_VIEWING_WIDTH) + (focalLocation.x - (WORLD_VIEWING_WIDTH / 2.0)) + float(i) * ALIASING_EPSILON/float(ALIASING_AMOUNT);
      				rayStartingLocation.y = (screenPosition.y * (WORLD_VIEWING_WIDTH * ASPECT_RATIO)) + (focalLocation.y - ((WORLD_VIEWING_WIDTH * ASPECT_RATIO) / 2.0)) + float(j) * ALIASING_EPSILON/float(ALIASING_AMOUNT);
      				rayStartingLocation.z = focalLocation.z + CAM_DISTANCE_TO_SCREEN;

      				rayStartingLocation.xz = rotationMat * (rayStartingLocation.xz-focalLocation.xz) + focalLocation.xz;
      				rayDirection = normalize(rayStartingLocation - focalLocation);

      				vec3 colorOfRay = getRayColor(rayStartingLocation, rayDirection);

      				totalLightRed += 1.0/float(ALIASING_AMOUNT*ALIASING_AMOUNT) * colorOfRay[0];
      				totalLightGreen += 1.0/float(ALIASING_AMOUNT*ALIASING_AMOUNT) * colorOfRay[1];
      				totalLightBlue += 1.0/float(ALIASING_AMOUNT*ALIASING_AMOUNT) * colorOfRay[2];
  	   		  }
  	  	  }     
        	gl_FragColor = vec4(totalLightRed, totalLightGreen, totalLightBlue, 1.0);
        }

    )).link();

    quad << vec2(-.5, -.5) << vec2(.5, -.5) << vec2(-.5, .5) << vec2(.5, .5);
    quad.upload();
    layout.create(shader, quad).attribute<float>("vertex", 2).check();
}

static vec3 eye = vec3(.01, 0.0, 0.0);
static vec3 axisX = vec3(1, 0, 0);
static vec3 axisY = vec3(0, 1, 0);
static vec3 axisZ = vec3(0, 0, 1);

static vec3 lightDirection = vec3(1, 0.5, 1);

static bool left = false;
static bool right = false;
static bool up = false;
static bool down = false;

static bool W = false;
static bool A = false;
static bool S = false;
static bool D = false;
static bool E = false;
static bool R = false;
static bool F = false;
static bool G = false;

static int enableMoon = 0;
static int enableSnow = 0;
static int enableRipples = 0;
static int enableClouds = 0;


static float angleX = 0;
static float angleY = 0;
static float angleOfRotation = 0.0;
static float maxVisibilityDueToFog = 40.0;
static float timeElapsed = 0.0;
const float MIN_FOG_AMOUNT = 60.0;
const float FOG_INCREMENT_AMOUNT = 1.0;

const float lightIncrementAmount = 0.01;

// Jump Variables
static bool isJumping = false;
static float currentJumpHeight = 0.0;
static float INITIAL_VELOCITY = 28.0;
static float GRAVITY_CONSTANT = -40.0;
static float timeAtBeginningOfJump = 0.0;

static void draw() {
  glClear(GL_COLOR_BUFFER_BIT);
  shader.use();

  shader.uniform("eye", eye);
  shader.uniform("axisX", axisX);
  shader.uniform("axisY", axisY);
  shader.uniform("axisZ", axisZ);
  shader.uniformFloat("angleOfRotation", angleOfRotation);
  shader.uniformFloat("timeElapsed", timeElapsed);
  shader.uniform("lightDirection", lightDirection);
  shader.uniformFloat("maxVisibilityDueToFog", maxVisibilityDueToFog);
  shader.uniformInt("enableMoon", enableMoon);
  shader.uniformInt("enableSnow", enableSnow);
  shader.uniformInt("enableRipples", enableRipples);
  shader.uniformInt("enableClouds", enableClouds);
  shader.uniformFloat("currentJumpHeight", currentJumpHeight);

  layout.draw(GL_TRIANGLE_STRIP);
  shader.unuse();
  glutSwapBuffers();
}

static void keyboard(unsigned char key, int x, int y) {
  switch (key) {
    case 27: 

  	printf("%s", gluErrorString(glGetError()));
    
    
    exit(0); break; // Escape
    case 'w': case 'W': W = true; break;
    case 'a': case 'A': A = true; break;
    case 's': case 'S': S = true; break;
    case 'd': case 'D': D = true; break;
    case 'e': case 'E': E = true; break;
    case 'r': case 'R': R = true; break;
    case 'f': case 'F': F = true; break;
    case 'g': case 'G': G = true; break;
  }
}

static void keyboardUp(unsigned char key, int x, int y) {
  switch (key) {
    case 'w': case 'W': W = false; break;
    case 'a': case 'A': A = false; break;
    case 's': case 'S': S = false; break;

    case 'i': case 'I': 
      if(enableSnow == 0)
      {
        enableSnow = 1;
      }
      else
      {
        enableSnow = 0;
      }
    break;

    case 'd': case 'D': D = false; break;
    case 'e': case 'E': E = false; break;
    case 'r': case 'R': R = false; break;
    case 'f': case 'F': F = false; break;
    case 'g': case 'G': G = false; break;

    case 'h': case 'H': 
      if(enableRipples == 0)
      {
        enableRipples = 1;
      }
      else
      {
        enableRipples = 0;
      }
    break;

    case 'm': case 'M': 
    	if(enableMoon == 0)
    	{
    		enableMoon = 1;
    	}
    	else
    	{
    		enableMoon = 0;
    	}
    break;
    case 'c': case 'C':
      if(enableClouds == 0)
      {
        enableClouds = 1;
      }
      else
      {
        enableClouds = 0;
      }
    break;

    case 32: 
      if(isJumping == false)
      {
        isJumping = true;
        timeAtBeginningOfJump = timeElapsed;
      }
    break;

  }
}

static void specialKey(int key, int x, int y) {
  switch (key) {
    case 27: exit(0); break; // Escape
    case GLUT_KEY_LEFT: left = true; break;
    case GLUT_KEY_RIGHT: right = true; break;
    case GLUT_KEY_UP: up = true; break;
    case GLUT_KEY_DOWN: down = true; break; 
  }
}

static void specialKeyUp(int key, int x, int y) {
  switch (key) {
    case GLUT_KEY_LEFT: left = false; break;
    case GLUT_KEY_RIGHT: right = false; break;
    case GLUT_KEY_UP: up = false; break;
    case GLUT_KEY_DOWN: down = false; break;

  }
}

static void centerCursor() {
  int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
  int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
  glutWarpPointer(screenWidth / 2, screenHeight / 2);
}

static void tick() {

  timeElapsed+=.007;
  if(left) {
    angleOfRotation += .1;
  }
  if(right) {
    angleOfRotation -= .1;
  }
  if(up) {
    //moving forward in the z direction of our coordinate frame
    eye.x += -sin(angleOfRotation) * .1;
    eye.z += cos(angleOfRotation) * .1;
  }
 if(down) {
  //moving backward in the z direction of our coordinate frame
    eye.x += -sin(angleOfRotation) * -.1;
    eye.z += cos(angleOfRotation) * -.1;
  }

  if(A)
  {
    lightDirection = vec3(lightDirection.x - lightIncrementAmount, lightDirection.y, lightDirection.z);
  }


  if(W)
  {
    lightDirection = vec3(lightDirection.x + lightIncrementAmount, lightDirection.y, lightDirection.z);
  }


  if(S)
  {
    lightDirection = vec3(lightDirection.x, lightDirection.y - lightIncrementAmount, lightDirection.z);
  }

  if(E)
  {
    lightDirection = vec3(lightDirection.x, lightDirection.y + lightIncrementAmount, lightDirection.z);
  }

  if(D)
  {
    lightDirection = vec3(lightDirection.x, lightDirection.y, lightDirection.z - lightIncrementAmount);
  }
  
  if(R)
  {
    lightDirection = vec3(lightDirection.x, lightDirection.y, lightDirection.z + lightIncrementAmount);
  }
  
  if(F)
  {
  	if(maxVisibilityDueToFog + FOG_INCREMENT_AMOUNT < MIN_FOG_AMOUNT)
  	{
  		maxVisibilityDueToFog += FOG_INCREMENT_AMOUNT;
  	}
  }
  
  if(G)
  {
  	if(maxVisibilityDueToFog - FOG_INCREMENT_AMOUNT > 0.0)
  	{
  		maxVisibilityDueToFog -= FOG_INCREMENT_AMOUNT;
  	}
  }

  if(isJumping)
  {
    float jumpTime = (timeElapsed - timeAtBeginningOfJump) * 5.0;

    currentJumpHeight = (0.5 * GRAVITY_CONSTANT * jumpTime * jumpTime) + (INITIAL_VELOCITY * jumpTime);

    if(currentJumpHeight < 0.0)
    {
      currentJumpHeight = 0.0;
      isJumping = false;
    }
  }

  glutPostRedisplay();
}

int main(int argc, char *argv[]) {
  printf("rotate: W/A/S/D\n");
  printf("move: up/down/left/right\n");
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutCreateWindow("Example");
  glutDisplayFunc(draw);
  glutSpecialFunc(specialKey);
  glutSpecialUpFunc(specialKeyUp);
  glutKeyboardFunc(keyboard);
  glutKeyboardUpFunc(keyboardUp);
  glutIdleFunc(tick);
  glutFullScreen();
  setup();
  glutMainLoop();
  return 0;
}
