;; REPEATED PARTICLE TRACKING
;; DEFINITIONS
(define ninjections 2)

;; RUN THE INITIAL INJECTION
(display (format #f "\n\n--------------------\nINITIAL INJECTION\n--------------------\n\n"))
(ti-menu-load-string (format #f "/report/dpm-sample initial-injection () outlet () () no"))

;; LOOP THE INJECTIONS
(DO ((x 1 (+ 1 x))) ((> x (- ninjections 1)))
	(ti-menu-load-string (format #f "/define/injections/list-particles file-injection"))
	(display (format #f "\n\n--------------------\nIteration ~a of ~a\n--------------------\n\n" (+ x 1) ninjections))
	(ti-menu-load-string (format #f "/report/dpm-sample file-injection () outlet () () no"))
)