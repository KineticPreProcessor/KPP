;;; kpp.el --- Major mode for KPP (Kinetic PreProcessor) files  -*- lexical-binding: t; -*-

;; Author:     Rolf Sander <rolf.sander@mpic.de>
;; Version:    1.2
;; Keywords:   languages, kpp, chemistry

;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2 of the License, or
;; (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;;; Commentary:
;;
;; This package provides a major mode for editing KPP (Kinetic PreProcessor)
;; files.  KPP is a software tool that assists in developing and solving
;; chemical kinetics problems.
;;
;; Installation:
;;   Copy kpp.el to a directory in `load-path' and add
;;     (require 'kpp)
;;   to your Emacs init file (~/.emacs or ~/.emacs.d/init.el).
;;
;; File extensions automatically associated with kpp-mode: .kpp, .eqn
;;
;; Note on .def and .spc extensions:
;;   These are intentionally NOT registered by default because they are used
;;   by many other file types (C preprocessor, Win32 module definitions,
;;   SPICE netlists, R packages, etc.).  Uncomment the relevant lines near
;;   the bottom of this file if you want them.
;;
;; Known issue:
;;   ":" inside KPP inline comments (between reaction products and the rate
;;   constant) may confuse font-lock.

;;; Code:

;; ---------------------------------------------------------------------------
;; Customization group
;; ---------------------------------------------------------------------------

(defgroup kpp nil
  "Major mode for editing KPP (Kinetic PreProcessor) files."
  :group  'languages
  :prefix "kpp-")

;; ---------------------------------------------------------------------------
;; User-configurable variables
;; ---------------------------------------------------------------------------

(defcustom kpp-comment-prefix "// "
  ;; Fix #10: renamed from `kpp-comment-region' to avoid collision with the
  ;; function of the same name.
  ;; Fix #13: changed from `defvar' with leading `*' to `defcustom'.
  "String inserted by \\[kpp-comment-region] at the start of each line."
  :type  'string
  :group 'kpp)

;; ---------------------------------------------------------------------------
;; Syntax table  (defined once at load time; shared by all kpp buffers)
;; ---------------------------------------------------------------------------

(defvar kpp-mode-syntax-table
  ;; Fix #11: was recreated (and re-allocated) inside `kpp-mode' on every
  ;; buffer visit.
  (let ((st (copy-syntax-table)))
    (modify-syntax-entry ?_ "w" st)  ; underscore is part of a word
    st)
  "Syntax table for `kpp-mode'.")

;; ---------------------------------------------------------------------------
;; Font-lock keyword specification
;; ---------------------------------------------------------------------------

(defconst kpp-font-lock-keywords
  ;; Fix #9: was a bare `setq' on an unbound symbol; `defconst' documents
  ;; intent and prevents inadvertent modification.
  (list
   ;; Reaction: "species = products : rate_constant ;"
   ;; Highlight the LHS + products (group 1) as a constant.
   '("^\\([^=\n]*=[^:\n]*\\):[^;\n]*;" 1 font-lock-constant-face)

   ;; Equation tag, e.g. <R001> or <j_NO2>
   ;; Fix #2: [A-z] spans ASCII 65-122 and accidentally includes the
   ;; characters [\]^_` (91-96).  Use [A-Za-z] instead.
   '("<[A-Za-z0-9_#]+>"  0 font-lock-variable-name-face t)

   ;; Curly-brace block comment: {comment text}
   '("{[^}\n]*}"          0 font-lock-comment-face t)

   ;; Fortran-90 / KPP inline comment: ! comment to end of line
   '("!.*"                0 font-lock-comment-face t)

   ;; {@...} – alternative LaTeX text
   ;; '("{@[^}]+}"           0 font-lock-doc-face t)

   ;; uncertainty of rate coefficient
   ;; '("{§[^}]*}"           0 font-lock-builtin-face t)

   ;; {$...} – alternative LaTeX text
   ;; '("{\\$[^}]+}"         0 font-lock-doc-face t)

   ;; {&...} – BibTeX reference
   ;; '("{&[^}]+}"           0 font-lock-function-name-face t)

   ;; {%...} – marker tag, e.g. {%TrG}
   ;; Fix #2: same [A-z] → [A-Za-z] correction as above.
   ;; '("{%[A-Za-z0-9#]+}"   0 font-lock-type-face t)

   ;; KPP sections, commands, and fragments
   ;; (Tables 3, 13, and 17 in the KPP manual / thesis)
   ;;
   ;; ORDERING RULE: in Emacs regex alternation (\|) the leftmost match
   ;; wins.  Any keyword that is a strict prefix of another must therefore
   ;; appear AFTER the longer form:
   ;;
   ;;   Pair                  Original order   Fixed order
   ;;   #CHECK / #CHECKALL    CHECKALL first   ✓ already correct
   ;;   #LOOKAT / #LOOKATALL  LOOKATALL first  ✓ already correct
   ;;   #TRANSPORT / …ALL     …ALL first       ✓ already correct
   ;;   #USE / #USES          USE first        ✗ BUG → fixed below (#Fix 3)
   ;;
   ;; \> (end-of-word boundary) prevents a short keyword from matching
   ;; as a prefix of an unknown longer token.
   (cons (concat
          "\\("
          "#ATOMS\\|#AUTOREDUCE"
          "\\|#CHECKALL\\|#CHECK"
          "\\|#DEFFIX\\|#DEFRAD\\|#DEFVAR"
          "\\|#DOUBLE\\|#DRIVER\\|#DUMMYINDEX"
          "\\|#ENDINLINE\\|#EQNTAGS\\|#EQUATIONS"
          "\\|#FAMILIES\\|#FUNCTION\\|#GRAPH"
          "\\|#HESSIAN"
          "\\|#INCLUDE\\|#INITVALUES\\|#INITIALIZE\\|#INLINE"
          "\\|#INTEGRATOR\\|#INTFILE"
          "\\|#JACOBIAN"
          "\\|#LANGUAGE"
          "\\|#LOOKATALL\\|#LOOKAT\\|#LUMP"
          "\\|#MEX\\|#MINVERSION\\|#MODEL\\|#MONITOR"
          "\\|#REORDER\\|#RUN"
          "\\|#SETFIX\\|#SETRAD\\|#SETVAR"
          "\\|#SPARSEDATA\\|#STOCHASTIC\\|#STOICMAT"
          "\\|#TRANSPORTALL\\|#TRANSPORT"
          "\\|#UPPERCASEF90"
          "\\|#USES\\|#USE"           ; Fix #3: #USES before #USE
          "\\|#WRITE_ATM\\|#WRITE_MAT\\|#WRITE_OPT\\|#WRITE_SPC"
          "\\|#XGRID\\|#YGRID\\|#ZGRID"
          "\\)\\>")
         'font-lock-keyword-face)

   ;; LaTeX note
   ;; '("//.*"               0 font-lock-string-face t)

   ;; C++-style line comment: // comment
   '("^//.*"              0 font-lock-comment-face t)

   ;; Fixme
   ;; '("qqq"                0 font-lock-warning-face t)       

   )
  "Font-lock keyword specification for `kpp-mode'.")

;; ---------------------------------------------------------------------------
;; Keymap
;; (defined before `define-derived-mode'; the macro detects and reuses it)
;; ---------------------------------------------------------------------------

(defvar kpp-mode-map
  (let ((map (make-sparse-keymap)))
    ;; C-c ; → comment / uncomment region
    (define-key map "\C-c;" #'kpp-comment-region)
    ;; Fix #14: original hard-coded 8 spaces and had no docstring.
    ;; Using `tab-width' respects the user/project preference and updates
    ;; automatically when `tab-width' is changed locally.
    (define-key map (kbd "TAB")
      (lambda ()
        "Insert `tab-width' spaces (never a literal TAB character)."
        (interactive)
        (insert-char ?\s tab-width)))
    map)
  "Keymap for `kpp-mode'.")

;; ---------------------------------------------------------------------------
;; Comment / uncomment a region   (adapted from wave-comment-region)
;; ---------------------------------------------------------------------------

(defun kpp-comment-region (beg-region end-region arg)
  ;; Fix #10: the variable formerly named `kpp-comment-region' clashed with
  ;; this function name.  The variable is now `kpp-comment-prefix'.
  "Comment every line in the region.
Inserts `kpp-comment-prefix' at the start of every line between
BEG-REGION and END-REGION.  With non-nil ARG, uncomment instead."
  (interactive "*r\nP")
  (let ((end-region-mark (make-marker))
        (save-point      (point-marker)))
    (set-marker end-region-mark end-region)
    (goto-char beg-region)
    (beginning-of-line)
    (if (not arg)
        ;; --- comment the region ---
        (progn
          (insert kpp-comment-prefix)
          (while (and (= (forward-line 1) 0)
                      (< (point) end-region-mark))
            (insert kpp-comment-prefix)))
      ;; --- uncomment the region ---
      (let ((com (regexp-quote kpp-comment-prefix)))
        (when (looking-at com)
          (delete-region (point) (match-end 0)))
        (while (and (= (forward-line 1) 0)
                    (< (point) end-region-mark))
          (when (looking-at com)
            (delete-region (point) (match-end 0))))))
    (goto-char save-point)
    (set-marker end-region-mark nil)
    (set-marker save-point nil)))

;; ---------------------------------------------------------------------------
;; Major mode definition
;; ---------------------------------------------------------------------------

;;;###autoload
(define-derived-mode kpp-mode fundamental-mode "KPP"
  ;; Fix #8:  `define-derived-mode' replaces the manual boilerplate that
  ;;          was in the original:
  ;;            • calls `kill-all-local-variables'         (Fix #4)
  ;;            • sets `major-mode' and `mode-name'
  ;;            • calls `use-local-map' with kpp-mode-map
  ;;            • calls `run-mode-hooks'
  ;;            • installs :syntax-table automatically
  "Major mode for editing KPP (Kinetic PreProcessor) files.

Turning on `kpp-mode' runs the hook `kpp-mode-hook'.

Known issue: `:' inside KPP inline comments (between reaction products
and the rate constant) may confuse font-lock.

\\{kpp-mode-map}"
  :syntax-table kpp-mode-syntax-table  ; Fix #11: reuse pre-built table
  ;; Fix #12: `setq-local' replaces the two-step
  ;;          `(make-local-variable ...) (setq ...)' pattern.
  (setq-local font-lock-defaults '((kpp-font-lock-keywords) t t))
  ;; KPP uses { } as block-comment delimiters.
  (setq-local comment-start      "{ ")
  (setq-local comment-end        " }")
  (setq-local comment-start-skip "{+\\s-*")   ; added: needed by comment cmds
  ;; Disable automatic line wrapping.
  (auto-fill-mode 0)
  ;; Fix #1 + #5: the original called `(turn-on-font-lock)' (deprecated)
  ;; and also registered `font-lock-mode' (a toggle!) on `kpp-mode-hook'.
  ;; With the hook in place, `turn-on-font-lock' ran first (enabling it),
  ;; then the hook ran `font-lock-mode' with no argument (toggling it OFF).
  ;; Net result: font-lock ended up disabled.  Fix: call `font-lock-mode'
  ;; with explicit argument 1 here; do NOT add it to `kpp-mode-hook'.
  (font-lock-mode 1))

;; ---------------------------------------------------------------------------
;; Auto-mode associations
;; ---------------------------------------------------------------------------

;; Fix #7: `add-to-list' is idempotent; the original `(setq … (cons …))'
;; pattern appended duplicate entries every time the file was reloaded.

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.kpp\\'" . kpp-mode))
;;;###autoload
(add-to-list 'auto-mode-alist '("\\.eqn\\'" . kpp-mode))

;; Fix #17: .def and .spc are intentionally omitted.  Those extensions are
;; claimed by many unrelated file types (C preprocessor, Win32 .def, SPICE,
;; R packages, …).  Uncomment the lines below ONLY if you use these
;; extensions exclusively for KPP files.
;;
;; (add-to-list 'auto-mode-alist '("\\.def\\'" . kpp-mode))
;; (add-to-list 'auto-mode-alist '("\\.spc\\'" . kpp-mode))

(provide 'kpp)

;;; kpp.el ends here
