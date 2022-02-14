;; REPEATED PARTICLE TRACKING
;; Inertial microfluidic devices often have several repeated segments
;; to modify fluid flow and particle positions within a channel.
;; Meshing and simulating the entire device would be computationally expensive
;; and applying periodic boundary conditions assumes the device is infinitely long.
;; 
;; To simulate particle flow in an inertial microfluidic device with finite repeated steps,
;; we first release particles at the inlet from an initial injection DPM definition.
;; Particles are then sampled at the outlet plane and re-injected in subsequent steps.


;; DEFINITIONS
(define ninjections 30)


;; RUN THE INITIAL INJECTION
(display (format #f "\n\n--------------------\nINITIAL INJECTION\n--------------------\n\n"))
(ti-menu-load-string (format #f "/report/dpm-sample initial-injection () outlet () () no"))


;; LOOP THE INJECTIONS
(DO ((x 1 (+ 1 x))) ((> x (- ninjections 1)))
	(ti-menu-load-string (format #f "/define/injections/list-particles file-injection"))
	(display (format #f "\n\n--------------------\nIteration ~a of ~a\n--------------------\n\n" (+ x 1) ninjections))
	(ti-menu-load-string (format #f "/report/dpm-sample file-injection () outlet () () no"))
)