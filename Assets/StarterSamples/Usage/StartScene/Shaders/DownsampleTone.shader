Shader "Unlit/DownsampleTone"
{
    Properties
    {
        _ToneImg ("Texture", 2D) = "white" {}
        _TileSz ("Tile Size", int) = 0
        _IgnoreZeros ("Ignore Zeros", int) = 0
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

      Pass { // some shaders require multiple passes
         GLSLPROGRAM // here begins the part in Unity's GLSL


         #extension GL_ARB_texture_rectangle : enable
         #extension GL_EXT_gpu_shader4 : enable



         uniform sampler2D _ToneImg;
         uniform highp int _TileSz;
         uniform highp int _IgnoreZeros;


         #ifdef VERTEX // here begins the vertex shader

         out vec4 texCoord;
         void main() // all vertex shaders define a main() function
         {
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
               // this line transforms the predefined attribute
               // gl_Vertex of type vec4 with the predefined
               // uniform gl_ModelViewProjectionMatrix of type mat4
               // and stores the result in the predefined output
               // variable gl_Position of type vec4.
            texCoord = gl_MultiTexCoord0;
         }

         #endif // here ends the definition of the vertex shader


         #ifdef FRAGMENT // here begins the fragment shader


         in vec4 texCoord;
         void main() // all fragment shaders define a main() function
         {
            vec2 im_sz = vec2(textureSize(_ToneImg,0));
            float x, y;
            vec3 sum = vec3(0.);
            float count = 0.;
            for( y=texCoord.y;
            y<min(texCoord.y+float(_TileSz), im_sz.y);
             y++ ) {
               for( x=texCoord.x; x<min(texCoord.x+float(_TileSz), im_sz.x); x++ ) {
                     // x - tone before tone mapping
                     // y - tone^2
                     // z - tone value after tone mapping
                     vec3 tone = texture(_ToneImg, vec2(float(x),float(y)) ).xyz;
                     if( !bool(_IgnoreZeros) || tone.x > -12.9 ) {
                        sum += tone;
                        count++;
                     }
               }
             }
            if( bool(_IgnoreZeros) ) {
               sum = sum/(count+1.0);
            } else {
               sum /= count;
            }

            // Compute and store std in the y-channel
            // std = sqrt(E(x^2) - E(x)^2)
            sum.y = sqrt(max(sum.y - sum.x*sum.x,0.0));


            gl_FragColor.rgb = sum;
            gl_FragColor.a = 1.;


         }


         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
      }

    }
}
