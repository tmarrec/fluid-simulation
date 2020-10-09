#include <QKeyEvent>
#include <QTimer>
#include <QObject>

class GlWidget;

class InputManager : public QObject
{
Q_OBJECT
public:
	InputManager(GlWidget* __glWidget);
	void keyPressEvent(QKeyEvent *event);
	void keyReleaseEvent(QKeyEvent *event);
	
	GlWidget* _glWidget;

private:
	bool _moveFront = false;
	bool _moveBack = false;
	bool _moveLeft = false;
	bool _moveRight = false;
	bool _moveUp = false;
	bool _moveDown = false;
	QTimer* _input_timer;
	void _process_inputs();
};
