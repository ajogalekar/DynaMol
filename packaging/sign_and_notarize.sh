#!/bin/bash
# DynaMol: Developer ID sign -> (optional) notarize -> staple -> verify
#
# Prereqs:
#   1. Developer ID Application cert installed  (security find-identity -v -p codesigning)
#   2. For notarization: creds stored in keychain profile "dynamol-notary":
#        xcrun notarytool store-credentials dynamol-notary \
#          --apple-id <APPLE_ID> --team-id B82PSPSJ64 --password <APP_SPECIFIC_PW>
#
# Usage (sign only):
#   IDENTITY="Developer ID Application: Ashutosh Jogalekar (B82PSPSJ64)" \
#   APP=/path/to/DynaMol.app  bash packaging/sign_and_notarize.sh
# Usage (sign + notarize + staple): also pass DMG=/path/to/DynaMol.dmg
set -uo pipefail

IDENTITY="${IDENTITY:?set IDENTITY (see: security find-identity -v -p codesigning)}"
APP="${APP:?set APP to the DynaMol.app to sign}"
ENT="${ENT:-$(cd "$(dirname "$0")" && pwd)/dynamol.entitlements}"
NOTARY_PROFILE="${NOTARY_PROFILE:-dynamol-notary}"
DMG="${DMG:-}"

echo ">> identity: $IDENTITY"
echo ">> app:      $APP"
echo ">> entitlements: $ENT"

# Strip xattrs that break codesign
xattr -cr "$APP" 2>/dev/null || true

# Sign every nested Mach-O. Executables get the entitlements (they run as processes,
# e.g. the bundled python backend); shared libraries/bundles get hardened runtime only.
echo ">> signing nested Mach-O inside-out (a few minutes; ~900 files)..."
libs=0; exes=0; fail=0
while IFS= read -r -d '' f; do
  ft=$(file -b "$f" 2>/dev/null)
  case "$ft" in
    *Mach-O*executable*)
      if codesign --force --timestamp --options runtime --entitlements "$ENT" --sign "$IDENTITY" "$f" >/dev/null 2>&1; then
        exes=$((exes+1)); else fail=$((fail+1)); echo "   FAIL(exe): $f"; fi ;;
    *Mach-O*)
      if codesign --force --timestamp --options runtime --sign "$IDENTITY" "$f" >/dev/null 2>&1; then
        libs=$((libs+1)); else fail=$((fail+1)); echo "   FAIL(lib): $f"; fi ;;
  esac
done < <(find "$APP" -type f -print0)
echo "   signed: libs=$libs exes=$exes  failures=$fail"

# Sign nested bundles (frameworks / helper .apps), deepest first
while IFS= read -r -d '' b; do
  [ "$b" = "$APP" ] && continue
  codesign --force --timestamp --options runtime --sign "$IDENTITY" "$b" >/dev/null 2>&1 \
    || echo "   FAIL(bundle): $b"
done < <(find "$APP" \( -name "*.framework" -o -name "*.app" \) -print0 2>/dev/null | \
         awk 'BEGIN{RS="\0";ORS="\0"}{print length($0)"\t"$0}' | sort -zrn | cut -zf2-)

# Sign the top-level app with entitlements + hardened runtime
echo ">> signing main app bundle..."
codesign --force --timestamp --options runtime --entitlements "$ENT" --sign "$IDENTITY" "$APP"

# Verify
echo ">> verifying signature..."
codesign --verify --deep --strict --verbose=2 "$APP" 2>&1 | tail -4
codesign -dvv "$APP" 2>&1 | grep -Ei "Authority=Developer ID Application|TeamId|flags" | head -3

if [ -n "$DMG" ]; then
  echo ">> signing DMG + notarize + staple..."
  codesign --force --timestamp --sign "$IDENTITY" "$DMG"
  xcrun notarytool submit "$DMG" --keychain-profile "$NOTARY_PROFILE" --wait
  xcrun stapler staple "$DMG"
  spctl -a -vvv -t open --context context:primary-signature "$DMG" 2>&1 | head -4
  echo ">> DONE — ticket stapled."
else
  echo ">> App signed. Re-run with DMG=<path> once the DMG is built + creds stored."
fi
