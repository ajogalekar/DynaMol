import AppKit
let directory = URL(fileURLWithPath: CommandLine.arguments[1], isDirectory: true)
try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
for (points, scale) in [(16,1),(16,2),(32,1),(32,2),(128,1),(128,2),(256,1),(256,2),(512,1),(512,2)] {
    let pixels = points * scale
    let image = NSImage(size: NSSize(width: pixels, height: pixels))
    image.lockFocus()
    let transform = AffineTransform(scale: CGFloat(pixels) / 1024)
    (transform as NSAffineTransform).concat()
    let background = NSBezierPath(roundedRect: NSRect(x: 60, y: 60, width: 904, height: 904), xRadius: 215, yRadius: 215)
    let gradient = NSGradient(starting: NSColor(srgbRed: 0.12, green: 0.30, blue: 0.28, alpha: 1), ending: NSColor(srgbRed: 0.035, green: 0.10, blue: 0.15, alpha: 1))!
    gradient.draw(in: background, angle: -60)
    NSColor(srgbRed: 0.51, green: 0.90, blue: 0.77, alpha: 1).setStroke()
    for degrees in [0.0, 60.0, 120.0] {
        let orbit = NSBezierPath(ovalIn: NSRect(x: 192, y: 370, width: 640, height: 284))
        var rotation = AffineTransform(translationByX: 512, byY: 512)
        rotation.rotate(byDegrees: degrees)
        rotation.translate(x: -512, y: -512)
        orbit.transform(using: rotation)
        orbit.lineWidth = 22
        orbit.stroke()
    }
    NSColor(srgbRed: 0.68, green: 0.96, blue: 0.86, alpha: 1).setFill()
    NSBezierPath(ovalIn: NSRect(x: 464, y: 464, width: 96, height: 96)).fill()
    image.unlockFocus()
    let bitmap = NSBitmapImageRep(data: image.tiffRepresentation!)!
    let png = bitmap.representation(using: .png, properties: [:])!
    let name = "icon_\(points)x\(points)\(scale == 2 ? "@2x" : "").png"
    try png.write(to: directory.appendingPathComponent(name))
}
