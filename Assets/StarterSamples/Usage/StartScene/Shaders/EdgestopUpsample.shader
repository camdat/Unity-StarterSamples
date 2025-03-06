Shader "Unlit/EdgestopUpsample"
{
    Properties
    {
        _ToneImgDs ("Texture", 2D) = "white" {}
        _Image ("Edge Stop Texture", 2D) = "white" {}
        _Sigma ("Sigma", float) = 1
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

      Pass { // some shaders require multiple passes
         GLSLPROGRAM // here begins the part in Unity's GLSL


         #extension GL_ARB_texture_rectangle : enable
         #extension GL_EXT_gpu_shader4 : enable



         uniform sampler2D _ToneImgDs; // Downsampled tone image

         uniform sampler2D _Image; // Guide texture used for edge-stopping

         uniform float _Sigma;  // Spread of the kernel


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
         vec4 range_kernel( vec4 d )
         {
            const float epsilon = 0.5;
            return max( vec4(epsilon), 1.-abs(d)/_Sigma );
         }

         void main()
         {
            // 2DRect textures have texels shifted by 0.5
            // vec2 tex_cord_rect = texCoord.st; // + vec2(-1,0);
            // vec4 im = vec4( texture(_ToneImgDs, tex_cord_rect - vec2(0,0)).z,
            //                texture(_ToneImgDs, tex_cord_rect - vec2(1,0)).z,
            //                texture(_ToneImgDs, tex_cord_rect - vec2(0,1)).z,
            //                texture(_ToneImgDs, tex_cord_rect - vec2(1,1)).z );

            // vec2 st = fract( tex_cord_rect.st );

            // // Spatial filter weights
            // vec4 weight_s = vec4( st.s*st.t, (1.-st.s)*st.t, st.s*(1.-st.t), (1.-st.s)*(1.-st.t) );

            // // Guidance texture value
            // float guide_v = texture(_Image, texCoord.xy).z; // Querry the tone value after tone mapping (channel "z")

            // // Range filter weights
            // vec4 weight_r = range_kernel( im-guide_v );

            // vec4 weight = weight_s*weight_r; // Total weights

            // float v = dot( im, weight ) / (weight.x+weight.y+weight.z+weight.w);
            gl_FragColor.rgb = texture(_Image, texCoord.st).rgb;
            // gl_FragColor.rgb = vec3( v );
            gl_FragColor.a = 1.;
         }



         #endif // here ends the definition of the fragment shader

         ENDGLSL // here ends the part in GLSL
      }

    }
}
